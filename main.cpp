#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>

#include <array>
#include <iostream>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <map>
#include <set>
#include <thread>
#include <mutex>
#include <queue>
#include <functional>

#define range(c) c.begin(), c.end()

using Clock = std::chrono::high_resolution_clock;

static constexpr int SIDE_LENGTH = 100;

static const auto PEACH_PUFF = sf::Color(0xFFDAB9FF);
static const auto DARK_GREY = sf::Color(0x6B6E76FF);
static const auto SCARLET = sf::Color(0xfc2847ff);

static std::mt19937 gen;
static std::uniform_int_distribution<> dis;
static auto COLORS = std::vector{sf::Color::Cyan, sf::Color::Red, sf::Color::Blue};

static constexpr std::array<std::pair<int, int>, 4> SIDE_NEIGHBORS{{
    {-1, 0},
    {0, 1},
    {1, 0},
    {0, -1},
}};

using Id = std::size_t;

struct Input {
    // Grid indices
    sf::Vector2i coordinates;

    bool              isLeft;
    bool              isSet = false;
};

struct Entity;
struct Visual;

class Game {
public:
    Game(sf::RenderWindow& window);

    void CreateRoad(sf::Vector2i coordinates);
    void CreateSource(sf::Vector2i gridCoordinates);
    void CreateTarget(sf::Vector2i gridCoordinates);

    void RemoveRoad(sf::Vector2i coordinates);

    std::vector<sf::Vector2i> FindPath(Id source, Id target) const;

    const std::vector<Id>& Get(sf::Vector2i gridCoordinates) const;

    std::shared_ptr<Entity> Get(Id id) const;

    std::vector<Id> GetFilledTargets() const;

    std::vector<sf::Vector2i> GetFreeCells() const;

    bool Has(sf::Vector2i gridCoordinates) const;

    bool IsOccupied(sf::Vector2i gridCoordinates) const;
    bool IsOccupied(int x, int y) const;

    bool IsRoad(sf::Vector2i coordinates) const;

    void ProcessInput();

    void HandleInput(Input input);

    void Render(double part);

    sf::Vector2i FromGrid(sf::Vector2i coordinates) const;
    sf::Vector2i FromGrid(int x, int y) const;
    sf::Vector2i ToGrid(sf::Vector2i coordinates) const;
    sf::Vector2i ToGrid(int x, int y) const;

    void Update();

    uint32_t GetHeight() const;

    uint32_t GetWidth() const;

    template <class T, class ...Args>
    std::shared_ptr<T> AddEntity(sf::Vector2i gridCoordinates, Args&& ...args);

    void RegisterToRender(Visual visual);

    void AddConnection(Id from, Id to);

    const std::unordered_set<Id>& RoadConnections(Id from);

private:
    sf::RenderWindow& window_;

    uint32_t cellByWidth_;
    uint32_t xDiff_;
    uint32_t cellByHeight_;
    uint32_t yDiff_;

    std::vector<std::vector<std::vector<Id>>> grid_;
    std::vector<std::shared_ptr<Entity>>      entities_;
    std::vector<Visual> toRender_;

    std::vector<Id> sources_;
    std::vector<Id> targets_;

    struct RepeatingTask {
        RepeatingTask(Clock::duration interval, std::function<void()> handler)
            : tp(Clock::now())
            , interval(interval)
            , handler(std::move(handler))
        {}

        Clock::time_point tp;
        Clock::duration interval;
        std::function<void()> handler;

        bool operator<(const RepeatingTask& rhs) const {
            return tp > rhs.tp;
        }
    };

    std::priority_queue<RepeatingTask> repeatingOperations_;

//    std::vector<std::thread> repeatingOperations_;
//    std::mutex mutex_;
//    bool                                      paused_ = false;
    Id id_ = 0;
    std::unordered_map<Id, std::unordered_set<Id>> roadGraph_;
};

struct Entity {
    Entity(Id id, sf::Vector2i position) : id(id), gridPosition(position) {}

    virtual ~Entity() = default;

    virtual void Update(Game& game) = 0;

    sf::Vector2i gridPosition;
    Id                id;
    bool removed = false;
};

struct Visual {
    Visual(std::unique_ptr<sf::Shape> shape) : shape(std::move(shape)) {}

    Visual(Visual&&) = default;
    Visual& operator=(Visual&&) = default;

    virtual ~Visual() = default;

    virtual void Render(sf::RenderWindow& window, double part) const {
        window.draw(*shape);
    }

    std::unique_ptr<sf::Shape> shape;
    sf::Color color = COLORS.at(dis(gen) % COLORS.size());
};

struct VisualEntity : public Entity, public Visual {
    VisualEntity(Id id, sf::Vector2i position, std::unique_ptr<sf::Shape> shape)
        : Entity(id, position)
        , Visual(std::move(shape))
    {}

    void SetPosition(float x, float y) {
        shape->setPosition({x, y});
    }
};

sf::CircleShape CAR_SHAPE{10};
sf::RectangleShape ROAD_SHAPE{sf::Vector2f{34, 34}};
sf::RectangleShape SOURCE_SHAPE{sf::Vector2f{50, 50}};
sf::RectangleShape TARGET_SHAPE{sf::Vector2f{90, 90}};

struct Car : public VisualEntity {
    Car(Id id, sf::Vector2i position, Id source)
        : VisualEntity(id, position, std::make_unique<sf::CircleShape>(CAR_SHAPE))
        , source(source)
    {}

    virtual void Update(Game& game) override {
        if (path.empty()) {
            auto targets = game.GetFilledTargets();
            if (targets.empty()) {
                return;
            }
            auto targetId = targets[dis(gen) % targets.size()];
            path = game.FindPath(source, targetId);
        }

        auto position = shape->getPosition();
        shape->move(speed);
    }

    static constexpr float FAST_SPEED = 10;
    static constexpr float MID_SPEED = 6;
    static constexpr float SLOW_SPEED = 2;
    sf::Vector2f speed;
    std::vector<sf::Vector2i> path;
    std::size_t curInd = 0;
    Id source;
};

struct V2iComp {
    bool operator()(sf::Vector2i lhs, sf::Vector2i rhs) const {
        return std::tie(lhs.x, lhs.y) < std::tie(rhs.x, rhs.y);
    }
};

struct Road : public VisualEntity {
    Road(Id id, sf::Vector2i position) : VisualEntity(id, position, std::make_unique<sf::RectangleShape>(ROAD_SHAPE)) {}

    void AddSegment(Game& game, sf::Vector2i direction) {
        if (!segments.contains(direction)) {
            auto segmentShape = std::make_unique<sf::RectangleShape>(
                static_cast<sf::RectangleShape&>(*shape));
            auto [width, height] = segmentShape->getSize();
            segmentShape->move(
                direction.x * (SIDE_LENGTH - width) / 2,
                direction.y * (SIDE_LENGTH - height) / 2);
            segments.try_emplace(direction, std::move(segmentShape));
        }
    }

    virtual void Render(sf::RenderWindow& window, double part) const {
        Visual::Render(window, part);
        for (const auto& [_, visual] : segments) {
            visual.Render(window, part);
        }
    }

    void Connect(Game& game) {
        for (auto [dx, dy] : SIDE_NEIGHBORS) {
            const sf::Vector2i coords = {gridPosition.x + dx, gridPosition.y + dy};
            if (game.Has(coords)) {
                const auto& neighbors = game.Get(coords);
                for (const auto& neighborId : neighbors) {
                    auto neighbor = game.Get(neighborId);
                    if (!neighbor->removed) {
                        if (auto road = std::dynamic_pointer_cast<Road>(neighbor)) {
                            if (!game.RoadConnections(id).contains(road->id)) {
                                game.AddConnection(id, road->id);
                                AddSegment(game, {dx, dy});
                                road->AddSegment(game, {-dx, -dy});
                            }
                        }
                    }
                }
            }
        }
    }

    virtual void Update(Game& game) override {
    }

    std::map<sf::Vector2i, Visual, V2iComp> segments;
};

struct Source : public VisualEntity {
    Source(Id id, sf::Vector2i position)
        : VisualEntity(id, position, std::make_unique<sf::RectangleShape>(SOURCE_SHAPE))
    {
        shape->setFillColor(SCARLET);
    }

    virtual void Update(Game& game) override {
    }
};

struct Target : public VisualEntity {
    Target(Id id, sf::Vector2i position)
        : VisualEntity(id, position, std::make_unique<sf::RectangleShape>(TARGET_SHAPE))
    {
        shape->setFillColor(SCARLET);
    }

    virtual void Update(Game& game) override {
    }

    std::size_t people = 0;
};

struct Config {
    uint32_t width             = 1280;
    uint32_t height            = 720;
    uint32_t antialiasingLevel = 16;
    uint32_t usPerUpdate       = 30'000;
};

Config LoadConfig(const std::filesystem::path& path) {
    std::ifstream file(path);
    Config        config;
    file >> config.width >> config.height >> config.antialiasingLevel >> config.usPerUpdate;
    return config;
}

void PlayGame(Config config, sf::RenderWindow& window) {
    Game      game(window);
    sf::Clock clock;
    sf::Int64 lag = 0;
    while (window.isOpen()) {
        window.clear(PEACH_PUFF);
        auto elapsed = clock.restart();
        lag += elapsed.asMicroseconds();

        game.ProcessInput();

        while (lag >= config.usPerUpdate) {
            game.Update();
            lag -= config.usPerUpdate;
        }

        game.Render(static_cast<double>(lag) / config.usPerUpdate);
        window.display();
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Config path is none" << std::endl;
        return 1;
    }

    auto config = LoadConfig(argv[1]);
    sf::ContextSettings settings;
    settings.antialiasingLevel = config.antialiasingLevel;
    sf::RenderWindow window({config.width, config.height}, "Game", sf::Style::Default, settings);

    window.setKeyRepeatEnabled(false);

    PlayGame(config, window);

    return 0;
}

Game::Game(sf::RenderWindow& window)
    : window_(window)
{
    entities_.resize(1024);
    cellByWidth_ = GetWidth() / SIDE_LENGTH;
    xDiff_ = GetWidth() % SIDE_LENGTH / 2;
    cellByHeight_ = GetHeight() / SIDE_LENGTH;
    yDiff_ = GetHeight() % SIDE_LENGTH / 2;

    for (std::size_t x = xDiff_; x < GetWidth(); x += SIDE_LENGTH) {
        auto shape = std::make_unique<sf::RectangleShape>(
            sf::Vector2f{1, static_cast<float>(cellByHeight_) * SIDE_LENGTH});
        shape->setPosition(x, yDiff_);
        shape->setFillColor(DARK_GREY);
        toRender_.emplace_back(std::move(shape));
    }

    for (std::size_t y = yDiff_; y < GetHeight(); y += SIDE_LENGTH) {
        auto shape = std::make_unique<sf::RectangleShape>(
            sf::Vector2f{static_cast<float>(cellByWidth_) * SIDE_LENGTH, 1});
        shape->setPosition(xDiff_, y);
        shape->setFillColor(DARK_GREY);
        toRender_.emplace_back(std::move(shape));
    }

    grid_.resize(cellByWidth_, {cellByHeight_, std::vector<Id>()});

    repeatingOperations_.emplace(std::chrono::seconds(4), [this]() {
        auto freeCells = GetFreeCells();
        if (freeCells.empty()) {
            return;
        }
        auto cell = freeCells[dis(gen) % freeCells.size()];
        CreateSource(cell);
    });

    repeatingOperations_.emplace(std::chrono::seconds(10), [this]() {
        auto freeCells = GetFreeCells();
        if (freeCells.empty()) {
            return;
        }
        auto cell = freeCells[dis(gen) % freeCells.size()];
        CreateTarget(cell);
    });

    repeatingOperations_.emplace(std::chrono::seconds(1), [this]() {
        if (targets_.empty()) {
            return;
        }
        auto targetId = targets_[dis(gen) % targets_.size()];
        ++static_cast<Target&>(*entities_.at(targetId)).people;
    });
}

void Game::ProcessInput() {
    Input     input;
    sf::Event event;
    if (window_.pollEvent(event)) {
        if (event.type == sf::Event::EventType::Closed) {
            std::exit(0);
        }
        if (event.type == sf::Event::KeyPressed) {
//            if (event.key.code == sf::Keyboard::Key::Escape) {
//                paused_ = !paused_;
//            }
        } else if (event.type == sf::Event::MouseButtonPressed) {
            if (event.mouseButton.button == sf::Mouse::Button::Left) {
                input.coordinates = {event.mouseButton.x, event.mouseButton.y};
                input.isLeft      = true;
                input.isSet       = true;
            } else if (event.mouseButton.button == sf::Mouse::Button::Right) {
                input.coordinates = {event.mouseButton.x, event.mouseButton.y};
                input.isLeft      = false;
                input.isSet       = true;
            }
        }
    }

    if (input.isSet) {
        HandleInput(input);
    }
}

void Game::HandleInput(Input input) {
    if (input.isLeft) {
        if (!IsOccupied(input.coordinates.x, input.coordinates.y)) {
            CreateRoad(input.coordinates);
        }
    } else if (IsRoad(input.coordinates)) {
        RemoveRoad(input.coordinates);
    }
}

void Game::Render(double part) {
    for (const auto& entity : entities_) {
        if (!entity || entity->removed) {
            continue;
        }
        if (auto ptr = std::dynamic_pointer_cast<Visual>(entity)) {
            ptr->Render(window_, part);
        }
    }
    for (const auto& vis : toRender_) {
        window_.draw(*vis.shape);
    }
}

void Game::Update() {
//    if (paused_) {
//        sf::sleep(sf::milliseconds(30));
//        return;
//    }
    auto top = repeatingOperations_.top();
    if (top.tp <= Clock::now()) {
        repeatingOperations_.pop();
        top.handler();
        top.tp += top.interval;
        repeatingOperations_.push(std::move(top));
    }

    auto size = entities_.size();
    for (std::size_t i = 0; i < size; ++i) {
//        auto createdObjects = entities_.at(i)->GetCreatedObjects();
//        entities_.insert(entities_.end(), createdObjects.begin(), createdObjects.end());
        auto& entity = entities_.at(i);
        if (entity && !entity->removed) {
            entity->Update(*this);
        }
    }
}

uint32_t Game::GetHeight() const {
    return window_.getSize().y;
}

uint32_t Game::GetWidth() const {
    return window_.getSize().x;
}

void Game::CreateRoad(sf::Vector2i coordinates) {
    if (IsRoad(coordinates)) {
        return;
    }
    auto [i, j] = ToGrid(coordinates);
    auto road = AddEntity<Road>({i, j}, coordinates);
    road->SetPosition(
        xDiff_ + i * SIDE_LENGTH + (SIDE_LENGTH - ROAD_SHAPE.getSize().x) / 2,
        yDiff_ + j * SIDE_LENGTH + (SIDE_LENGTH - ROAD_SHAPE.getSize().y) / 2);
    road->Connect(*this);
}

void Game::CreateSource(sf::Vector2i gridCoordinates) {
    auto position = FromGrid(gridCoordinates);
    auto source = AddEntity<Source>({gridCoordinates.x, gridCoordinates.y}, position);
    source->SetPosition(
        position.x + (SIDE_LENGTH - SOURCE_SHAPE.getSize().x) / 2,
        position.y + (SIDE_LENGTH - SOURCE_SHAPE.getSize().y) / 2);
    sources_.push_back(source->id);
}

void Game::CreateTarget(sf::Vector2i gridCoordinates) {
    auto position = FromGrid(gridCoordinates);
    auto target = AddEntity<Target>({gridCoordinates.x, gridCoordinates.y}, position);
    target->SetPosition(
        position.x + (SIDE_LENGTH - TARGET_SHAPE.getSize().x) / 2,
        position.y + (SIDE_LENGTH - TARGET_SHAPE.getSize().y) / 2);
    targets_.push_back(target->id);
}

void Game::RemoveRoad(sf::Vector2i coordinates) {
    const auto gridCoords = ToGrid(coordinates);
    for (auto id : Get(gridCoords)) {
        if (auto ptr = std::dynamic_pointer_cast<Road>(entities_.at(id))) {
            if (!ptr->removed) {
                ptr->removed = true;
                return;
            }
        }
    }
}

bool Game::IsRoad(sf::Vector2i coordinates) const {
    const auto gridCoords = ToGrid(coordinates);
    for (auto id : Get(gridCoords)) {
        if (auto ptr = std::dynamic_pointer_cast<Road>(entities_.at(id))) {
            if (!ptr->removed) {
                return true;
            }
        }
    }
    return false;
}

sf::Vector2i Game::FromGrid(sf::Vector2i coordinates) const {
    return FromGrid(coordinates.x, coordinates.y);
}

sf::Vector2i Game::FromGrid(int x, int y) const {
    return sf::Vector2i(
        xDiff_ + x * SIDE_LENGTH,
        yDiff_ + y * SIDE_LENGTH);
}

sf::Vector2i Game::ToGrid(sf::Vector2i coordinates) const {
    return ToGrid(coordinates.x, coordinates.y);
}

sf::Vector2i Game::ToGrid(int x, int y) const {
    return sf::Vector2i(
        x * grid_.size() / GetWidth(),
        y * grid_.at(0).size() / GetHeight());
}

template <class T, class... Args>
std::shared_ptr<T> Game::AddEntity(sf::Vector2i gridCoordinates, Args&& ... args) {
    auto entity = std::make_shared<T>(id_, std::forward<Args>(args)...);
    entity->gridPosition = gridCoordinates;

    entities_.at(id_) = entity;
    grid_.at(gridCoordinates.x).at(gridCoordinates.y).push_back(id_);
    ++id_;
    return entity;
}

void Game::RegisterToRender(Visual visual) {
    toRender_.emplace_back(std::move(visual));
}

const std::vector<Id>& Game::Get(sf::Vector2i gridCoordinates) const {
    const auto [x, y] = gridCoordinates;
    return grid_.at(x).at(y);
}

std::vector<sf::Vector2i> Game::GetFreeCells() const {
    std::vector<sf::Vector2i> freeCells;
    for (int i = 0; i < grid_.size(); ++i) {
        for (int j = 0; j < grid_[0].size(); ++j) {
            if (!IsOccupied(sf::Vector2i{i, j})) {
                freeCells.emplace_back(i, j);
            }
        }
    }
    return freeCells;
}

std::shared_ptr<Entity> Game::Get(Id id) const {
    return entities_.at(id);
}

bool Game::IsOccupied(sf::Vector2i gridCoordinates) const {
    const auto& entities = Get(gridCoordinates);
    return std::any_of(range(entities), [&](const auto& id) {
        return !entities_.at(id)->removed;
    });
}

bool Game::IsOccupied(int x, int y) const {
    return IsOccupied(ToGrid(x, y));
}

bool Game::Has(sf::Vector2i gridCoordinates) const {
    const auto [x, y] = gridCoordinates;
    return 0 <= x && x < grid_.size() && 0 <= y && y < grid_[0].size();
}

std::vector<Id> Game::GetFilledTargets() const {
    std::vector<Id> filled;
    for (auto targetId : targets_) {
        const auto& target = static_cast<Target&>(*entities_.at(targetId));
        if (target.people > 0) {
            filled.push_back(targetId);
        }
    }
    return filled;
}

std::vector<sf::Vector2i> Game::FindPath(Id source, Id target) const {
    std::vector<sf::Vector2i> path;

    using DistId = std::pair<std::size_t, std::size_t>;
    std::priority_queue<DistId, std::vector<DistId>, std::greater<>> queue;
    std::unordered_set<Id> visited;
    queue.emplace(0, source);
    visited.insert(source);
    while (!queue.empty()) {
        auto [curDist, id] = queue.top();
        auto [x, y] = entities_.at(id)->gridPosition;
        for (auto [dx, dy] : SIDE_NEIGHBORS) {
            sf::Vector2i neighborGridCoordinates{x + dx, y + dy};
            for (auto roadId : Get(neighborGridCoordinates)) {
                if (auto ptr = std::dynamic_pointer_cast<Road>(entities_.at(roadId))) {

                }
            }
            queue.emplace();
        }
    }

    return path;
}

void Game::AddConnection(Id from, Id to) {
    roadGraph_[from].insert(to);
    roadGraph_[to].insert(from);
}

const std::unordered_set<Id>& Game::RoadConnections(Id from) {
    return roadGraph_[from];
}
