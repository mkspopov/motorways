#pragma clang diagnostic push
#pragma ide diagnostic ignored "cppcoreguidelines-pro-type-static-cast-downcast"
#pragma ide diagnostic ignored "cert-msc51-cpp"
#pragma ide diagnostic ignored "cert-err58-cpp"
#pragma ide diagnostic ignored "modernize-use-override"

#include <SFML/Graphics.hpp>
#include <SFML/Window.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <mutex>
#include <queue>
#include <random>
#include <ranges>
#include <set>
#include <unordered_map>
#include <unordered_set>

////////////////////////////////////////////////////////////
template <typename T>
inline sf::Vector2<T> operator +(const sf::Vector2<T>& left, T right)
{
    return sf::Vector2<T>(left.x + right, left.y + right);
}

////////////////////////////////////////////////////////////
template <typename T>
inline sf::Vector2<T> operator +(T left, const sf::Vector2<T>& right)
{
    return right + left;
}

////////////////////////////////////////////////////////////
template <typename T>
inline sf::Vector2<T> operator -(const sf::Vector2<T>& left, T right)
{
    return sf::Vector2<T>(left.x - right, left.y - right);
}

////////////////////////////////////////////////////////////
template <typename T>
inline sf::Vector2<T> operator -(T left, const sf::Vector2<T>& right)
{
    return -(right - left);
}

////////////////////////////////////////////////////////////
template <typename T>
inline std::ostream& operator <<(std::ostream& os, const sf::Vector2<T>& right)
{
    return os << right.x << ' ' << right.y;
}

#define range(c) c.begin(), c.end()

using Clock = std::chrono::high_resolution_clock;

struct V2iComp {
    bool operator()(sf::Vector2i lhs, sf::Vector2i rhs) const {
        return std::tie(lhs.x, lhs.y) < std::tie(rhs.x, rhs.y);
    }
};

static constexpr int SIDE_LENGTH = 100;

static const auto PEACH_PUFF = sf::Color(0xFFDAB9FF);
static const auto DARK_GREY = sf::Color(0x6B6E76FF);
static const auto SCARLET = sf::Color(0xfc2847ff);

static constexpr auto SEED = 0;
static std::mt19937 gen(SEED);
static std::uniform_int_distribution<> dis;
static std::vector<sf::Color> COLORS;

sf::Font PURISA_FONT;

using Id = std::size_t;

// Different types for coordinates for compile-time error checking.
using GridXy = sf::Vector2i;
using WindXy = sf::Vector2f;

static const inline std::array<GridXy, 4> SIDE_NEIGHBORS{{
    {-1, 0},
    {0, 1},
    {1, 0},
    {0, -1},
}};

struct Input {
    WindXy coordinates;
    bool   isLeft;
    bool   isSet = false;
};

struct Entity;
struct Visual;

namespace render {

enum Priority {
    Lowest,
    Low,
    High,
    Size,
};

}

class Game {
public:
    Game(sf::RenderWindow& window);

    void CreateRoad(WindXy coordinates);
    void CreateSource(GridXy coordinates);
    void CreateTarget(GridXy coordinates);

    void RemoveRoad(GridXy coordinates);

    const std::vector<Id>& Get(GridXy coordinates) const;

    std::shared_ptr<Entity> Get(Id id) const;

    std::vector<Id> GetFilledTargets() const;

    std::vector<GridXy> GetFreeCells() const;

    bool Has(GridXy coordinates) const;

    bool IsOccupied(GridXy coordinates) const;
    bool IsOccupied(int x, int y) const;

    template <class U>
    bool Is(GridXy coordinates) const;

    void ProcessInput();

    void HandleInput(Input input);

    void Render(float part);

    WindXy FromGrid(GridXy coordinates) const;
    WindXy FromGrid(int x, int y) const;
    GridXy ToGrid(WindXy coordinates) const;
    GridXy ToGrid(int x, int y) const;

    void Update();

    uint32_t GetHeight() const;

    uint32_t GetWidth() const;

    template <class T, class ...Args>
    std::shared_ptr<T> AddEntity(GridXy coordinates, Args&& ...args);

    template <class ...Args>
    void RegisterToRender(int priority, Args&&... args);

    void AddConnection(Id from, Id to);
    void RemoveConnection(Id from, Id to);
    const std::unordered_set<Id>& RoadConnections(Id from);

    const auto& GetTargets() const {
        return targets_;
    }

    bool IsIn(WindXy windXy, GridXy gridXy) const {
        return ToGrid(windXy) == gridXy;
    }

private:
    sf::RenderWindow& window_;

    uint32_t cellByWidth_;
    uint32_t xDiff_;
    uint32_t cellByHeight_;
    uint32_t yDiff_;

    static const inline std::vector<Id> EMPTY_ARRAY;
    std::vector<std::vector<std::vector<Id>>> grid_;
    std::vector<std::shared_ptr<Entity>> entities_;
    std::array<std::vector<std::shared_ptr<Visual>>, static_cast<int>(render::Priority::Size)> toRender_;

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

//    bool                                      paused_ = false;
    Id id_ = 0;
    std::unordered_map<Id, std::unordered_set<Id>> roadGraph_;
};

struct Entity {
    Entity(Id id, GridXy position) : id(id), gridPosition(position) {}

    virtual ~Entity() = default;

    virtual void Update(Game& game) = 0;

    GridXy gridPosition;
    Id                id;
    bool alive = true;
};

struct Visual {
    Visual(std::unique_ptr<sf::Shape> shape) : shape(std::move(shape)) {}

    Visual(Visual&&) = default;
    Visual& operator=(Visual&&) = default;

    virtual ~Visual() = default;

    virtual void Render(sf::RenderWindow& window, float part) const {
        window.draw(*shape);
    }

    std::unique_ptr<sf::Shape> shape;
    sf::Color color = COLORS.at(dis(gen) % COLORS.size());
};

struct VisualEntity : public Entity, public Visual {
    VisualEntity(Id id, GridXy position, std::unique_ptr<sf::Shape> shape)
        : Entity(id, position)
        , Visual(std::move(shape))
    {}

    render::Priority priority = render::Priority::Lowest;
};

sf::CircleShape CAR_SHAPE{10};
sf::RectangleShape ROAD_SHAPE{sf::Vector2f{34, 34}};
sf::RectangleShape SOURCE_SHAPE{sf::Vector2f{50, 50}};
sf::RectangleShape TARGET_SHAPE{sf::Vector2f{90, 90}};

WindXy GetCenter(const std::unique_ptr<sf::Shape>& shape) {
    auto pos = shape->getPosition();
    if (auto circ = dynamic_cast<const sf::CircleShape*>(shape.get())) {
        return pos + circ->getRadius();
    } else if (auto rect = dynamic_cast<const sf::RectangleShape*>(shape.get())) {
        return pos + rect->getSize() / 2.0f;
    }
    throw std::runtime_error("Unknown shape in GetCenter()");
}

template <class T>
T Length(sf::Vector2<T> v) {
    return std::sqrt(v.x * v.x + v.y * v.y);
}

bool Contains(const std::unique_ptr<sf::Shape>& shape, WindXy point) {
    if (auto circ = dynamic_cast<const sf::CircleShape*>(shape.get())) {
        auto center = GetCenter(shape);
        return Length(point - center) <= circ->getRadius();
    } else if (auto rect = dynamic_cast<const sf::RectangleShape*>(shape.get())) {
        return sf::Rect(shape->getPosition(), rect->getSize()).contains(point);
    }
    throw std::runtime_error("Unknown shape in GetCenter()");
}

struct Target : public VisualEntity {
    Target(Id id, GridXy position)
        : VisualEntity(id, position, std::make_unique<sf::RectangleShape>(TARGET_SHAPE))
    {
        priority = render::Priority::Low;
        shape->setFillColor(SCARLET);
    }

    virtual void Update(Game& game) override {
    }

    virtual void Render(sf::RenderWindow& window, float part) const override {
        Visual::Render(window, part);
        sf::Text text(std::to_string(people) + ", " + std::to_string(takenPeople), PURISA_FONT);
        text.setPosition(shape->getPosition());
        text.setCharacterSize(30);
        text.setStyle(sf::Text::Bold);
        text.setFillColor(sf::Color::Black);
        window.draw(text);
    }

    bool HasPeopleToTake() const {
        return people > takenPeople;
    }

    std::size_t people = 0;
    std::size_t takenPeople = 0;
};

struct Road : public VisualEntity {
    Road(Id id, GridXy position)
        : VisualEntity(id, position, std::make_unique<sf::RectangleShape>(ROAD_SHAPE))
    {
        priority = render::Priority::Low;
    }

    void AddSegment(Game& game, GridXy direction) {
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

    virtual void Render(sf::RenderWindow& window, float part) const override {
        Visual::Render(window, part);
        for (const auto& [_, visual] : segments) {
            visual.Render(window, part);
        }
    }

    void Connect(Game& game);

    virtual void Update(Game& game) override {
    }

    std::map<GridXy, Visual, V2iComp> segments;
};

struct Car : public VisualEntity {
    Car(Id id, GridXy position, Id source)
        : VisualEntity(id, position, std::make_unique<sf::CircleShape>(CAR_SHAPE))
        , source(source)
        , target(-1)
    {
        shape->setFillColor(sf::Color::Black);
        priority = render::Priority::High;
    }

    void Go(auto carPath, Id targetId);

    bool IsFree() const {
        return state_ == nullptr;
    }

    void Move(std::vector<GridXy>& path, std::size_t& pathInd, bool& entering_, Game& game);

    virtual void Render(sf::RenderWindow& window, float part) const override {
        auto sh = static_cast<const sf::CircleShape&>(*shape);
        sh.move(speed * part);
        window.draw(sh);
    }

    virtual void Update(Game& game) override;

    static constexpr float FAST_SPEED = 10;
    static constexpr float MID_SPEED = 6;
    static constexpr float SLOW_SPEED = 2;
    sf::Vector2f speed;
    Id source;
    Id target;

    struct S;
    S* state_ = nullptr;
};

struct Source : public VisualEntity {
    Source(Id id, GridXy position)
        : VisualEntity(id, position, std::make_unique<sf::RectangleShape>(SOURCE_SHAPE))
    {
        shape->setFillColor(SCARLET);
        priority = render::Priority::Low;
    }

    virtual void Update(Game& game) override {
        for (auto& car : cars_) {
            if (!car) {
                car = game.AddEntity<Car>(gridPosition, id);
                car->shape->setPosition(
                    GetCenter(shape) - static_cast<const sf::CircleShape*>(car->shape.get())->getRadius());
            }
        }
        if (!cars_[0]->IsFree() && !cars_[1]->IsFree()) {
            return;
        }
        for (auto targetId : game.GetTargets()) {
            auto& target = static_cast<Target&>(*game.Get(targetId));
            if (target.HasPeopleToTake()) {
                using DistId = std::pair<std::size_t, Id>;
                std::priority_queue<DistId, std::vector<DistId>, std::greater<>> queue;
                std::unordered_map<Id, std::size_t> distances;
                queue.emplace(0, id);
                distances.emplace(id, 0);
                std::unordered_map<Id, Id> parents;
                while (!queue.empty()) {
                    auto [curDist, curId] = queue.top();
                    queue.pop();
                    if (curDist > distances[curId]) {
                        continue;
                    }
                    for (auto entityId : game.RoadConnections(curId)) {
                        if (entityId == targetId) {
                            distances[entityId] = curDist + 1;
                            parents[entityId] = curId;
                            queue = {};
                            break;
                        }
                        auto road = std::dynamic_pointer_cast<Road>(game.Get(entityId));
                        if (road && road->alive) {
                            auto dist = curDist + 1;
                            if (!distances.contains(entityId) || distances[entityId] > dist) {
                                distances[entityId] = dist;
                                queue.emplace(dist, entityId);
                                parents[entityId] = curId;
                            }
                        }
                    }
                }

                if (distances.contains(targetId)) {
                    std::vector<GridXy> path;
                    auto curId = target.id;
                    while (curId != id) {
                        path.emplace_back(game.Get(curId)->gridPosition);
                        curId = parents[curId];
                    }

                    path.emplace_back(gridPosition);
                    std::reverse(range(path));

                    ++target.takenPeople;

                    std::size_t carIndex = 0;
                    if (!cars_[carIndex]->IsFree()) {
                        carIndex = 1;
                    }

                    cars_[carIndex]->Go(std::move(path), targetId);

                    if (!cars_[1 - carIndex]->IsFree()) {
                        return;
                    }
                }
            }
        }
    }

    std::array<std::shared_ptr<Car>, 2> cars_;
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

    PURISA_FONT.loadFromFile("/usr/share/fonts/opentype/tlwg/Purisa-Bold.otf");

    while (window.isOpen()) {
        window.clear(PEACH_PUFF);
        auto elapsed = clock.restart();
        lag += elapsed.asMicroseconds();

        game.ProcessInput();

        while (lag >= config.usPerUpdate) {
            game.Update();
            lag -= config.usPerUpdate;
        }

        game.Render(static_cast<float>(lag) / config.usPerUpdate);
        window.display();
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Config path is none" << std::endl;
        return 1;
    }

    COLORS = {sf::Color::Cyan, sf::Color::Red, sf::Color::Blue};

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
        RegisterToRender(render::Priority::Low, std::move(shape));
    }

    for (std::size_t y = yDiff_; y < GetHeight(); y += SIDE_LENGTH) {
        auto shape = std::make_unique<sf::RectangleShape>(
            sf::Vector2f{static_cast<float>(cellByWidth_) * SIDE_LENGTH, 1});
        shape->setPosition(xDiff_, y);
        shape->setFillColor(DARK_GREY);
        RegisterToRender(render::Priority::Low, std::move(shape));
    }

    grid_.resize(cellByWidth_, {cellByHeight_, std::vector<Id>()});

    repeatingOperations_.emplace(std::chrono::seconds(4), [this]() {
        auto freeCells = GetFreeCells();
        std::erase_if(freeCells, [this](auto gridXy) {
            for (auto dNeighbor : SIDE_NEIGHBORS) {
                if(Is<Target>(gridXy + dNeighbor)) {
                    return true;
                }
            }
            return false;
        });
        if (freeCells.empty()) {
            return;
        }
        auto cell = freeCells[dis(gen) % freeCells.size()];
        CreateSource(cell);
    });

    repeatingOperations_.emplace(std::chrono::seconds(10), [this]() {
        auto freeCells = GetFreeCells();
        std::erase_if(freeCells, [this](auto gridXy) {
            for (auto dNeighbor : SIDE_NEIGHBORS) {
                if(Is<Source>(gridXy + dNeighbor)) {
                    return true;
                }
            }
            return false;
        });
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
    sf::Event event; // NOLINT(cppcoreguidelines-pro-type-member-init)
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
                input.coordinates = WindXy(event.mouseButton.x, event.mouseButton.y);
                input.isLeft      = true;
                input.isSet       = true;
            } else if (event.mouseButton.button == sf::Mouse::Button::Right) {
                input.coordinates = WindXy(event.mouseButton.x, event.mouseButton.y);
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
    } else if (Is<Road>(ToGrid(input.coordinates))) {
        RemoveRoad(ToGrid(input.coordinates));
    }
}

void Game::Render(float part) {
    for (const auto& visuals : toRender_) {
        for (const auto& vis : visuals) {
            auto ptr = std::dynamic_pointer_cast<VisualEntity>(vis);
            if (!ptr || ptr->alive) {
                vis->Render(window_, part);
            }
        }
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
        auto& entity = entities_.at(i);
        if (entity && entity->alive) {
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

void Game::CreateRoad(WindXy coordinates) {
    const auto gridXy = ToGrid(coordinates);
    if (Is<Road>(gridXy)) {
        return;
    }
    auto road = AddEntity<Road>(gridXy);
    const auto windXy = FromGrid(gridXy);
    road->shape->setPosition(
        windXy.x + (SIDE_LENGTH - ROAD_SHAPE.getSize().x) / 2,
        windXy.y + (SIDE_LENGTH - ROAD_SHAPE.getSize().y) / 2);
    road->Connect(*this);
}

void Game::CreateSource(GridXy coordinates) {
    auto position = FromGrid(coordinates);
    auto source = AddEntity<Source>(coordinates);
    source->shape->setPosition(
        position.x + (SIDE_LENGTH - SOURCE_SHAPE.getSize().x) / 2,
        position.y + (SIDE_LENGTH - SOURCE_SHAPE.getSize().y) / 2);
    sources_.push_back(source->id);

    for (auto diff : SIDE_NEIGHBORS) {
        for (auto entityId : Get(coordinates + diff)) {
            if (auto road = std::dynamic_pointer_cast<Road>(Get(entityId))) {
                road->AddSegment(*this, -diff);
                AddConnection(source->id, road->id);
            }
        }
    }
}

void Game::CreateTarget(GridXy coordinates) {
    auto position = FromGrid(coordinates);
    auto target = AddEntity<Target>(coordinates);
    target->shape->setPosition(
        position.x + (SIDE_LENGTH - TARGET_SHAPE.getSize().x) / 2,
        position.y + (SIDE_LENGTH - TARGET_SHAPE.getSize().y) / 2);
    targets_.push_back(target->id);

    for (auto diff : SIDE_NEIGHBORS) {
        for (auto entityId : Get(coordinates + diff)) {
            if (auto road = std::dynamic_pointer_cast<Road>(Get(entityId))) {
                road->AddSegment(*this, -diff);
                AddConnection(target->id, road->id);
            }
        }
    }
}

void Game::RemoveRoad(GridXy coordinates) {
    for (auto id : Get(coordinates)) {
        if (auto ptr = std::dynamic_pointer_cast<Road>(entities_.at(id))) {
            if (ptr->alive) {
                ptr->alive = false;
                for (auto it = roadGraph_[id].begin(); it != roadGraph_[id].end(); ) {
                    roadGraph_[*it].erase(id);
                    it = roadGraph_[id].erase(it);
                }
                return;
            }
        }
    }
}

template <class U>
bool Game::Is(GridXy coordinates) const {
    return std::ranges::any_of(Get(coordinates), [this](auto id) {
        auto ptr = std::dynamic_pointer_cast<U>(entities_.at(id));
        return ptr && ptr->alive;
    });
}

WindXy Game::FromGrid(GridXy coordinates) const {
    return FromGrid(coordinates.x, coordinates.y);
}

WindXy Game::FromGrid(int x, int y) const {
    return {
        static_cast<float>(xDiff_ + x * SIDE_LENGTH),
        static_cast<float>(yDiff_ + y * SIDE_LENGTH)};
}

GridXy Game::ToGrid(WindXy coordinates) const {
    return ToGrid(coordinates.x, coordinates.y);
}

GridXy Game::ToGrid(int x, int y) const {
    return {
        static_cast<int>(x * grid_.size() / GetWidth()),
        static_cast<int>(y * grid_.at(0).size() / GetHeight())};
}

template <class T, class... Args>
std::shared_ptr<T> Game::AddEntity(GridXy coordinates, Args&& ... args) {
    auto entity = std::make_shared<T>(id_, coordinates, std::forward<Args>(args)...);
    entity->gridPosition = coordinates;

    entities_.at(id_) = entity;
    if (auto ptr = std::dynamic_pointer_cast<VisualEntity>(entity)) {
        RegisterToRender(ptr->priority, std::move(ptr));
    }
    grid_.at(coordinates.x).at(coordinates.y).push_back(id_);
    ++id_;
    return entity;
}

template <class ...Args>
void Game::RegisterToRender(int priority, Args&&... args) {
    if constexpr (std::is_same_v<Args..., std::shared_ptr<Visual>> ||
        std::is_same_v<Args..., std::shared_ptr<VisualEntity>>)
    {
        toRender_.at(priority).emplace_back(std::move(args...));
    } else {
        toRender_.at(priority).emplace_back(std::make_shared<Visual>(std::forward<Args>(args)...));
    }
}

const std::vector<Id>& Game::Get(GridXy coordinates) const {
    const auto [x, y] = coordinates;
    if (0 <= x && x < grid_.size() && 0 <= y && y < grid_[0].size()) {
        return grid_[x][y];
    }
    return EMPTY_ARRAY;
}

std::vector<GridXy> Game::GetFreeCells() const {
    std::vector<GridXy> freeCells;
    for (int i = 0; i < grid_.size(); ++i) {
        for (int j = 0; j < grid_[0].size(); ++j) {
            if (!IsOccupied(GridXy{i, j})) {
                freeCells.emplace_back(i, j);
            }
        }
    }
    return freeCells;
}

std::shared_ptr<Entity> Game::Get(Id id) const {
    return entities_.at(id);
}

bool Game::IsOccupied(GridXy coordinates) const {
    const auto& entities = Get(coordinates);
    return std::any_of(range(entities), [&](const auto& id) {
        return entities_.at(id)->alive;
    });
}

bool Game::IsOccupied(int x, int y) const {
    return IsOccupied(ToGrid(x, y));
}

bool Game::Has(GridXy coordinates) const {
    const auto [x, y] = coordinates;
    return 0 <= x && x < grid_.size() && 0 <= y && y < grid_[0].size();
}

std::vector<Id> Game::GetFilledTargets() const {
    std::vector<Id> filled;
    for (auto targetId : targets_) {
        const auto& target = static_cast<Target&>(*entities_.at(targetId));
        if (target.HasPeopleToTake()) {
            filled.push_back(targetId);
        }
    }
    return filled;
}

void Game::AddConnection(Id from, Id to) {
    roadGraph_[from].insert(to);
    roadGraph_[to].insert(from);
}

void Game::RemoveConnection(Id from, Id to) {
    roadGraph_[from].erase(to);
    roadGraph_[to].erase(from);
}

const std::unordered_set<Id>& Game::RoadConnections(Id from) {
    return roadGraph_[from];
}

void Road::Connect(Game& game) {
    for (auto [dx, dy] : SIDE_NEIGHBORS) {
        const GridXy coords = {gridPosition.x + dx, gridPosition.y + dy};
        const auto& neighbors = game.Get(coords);
        for (const auto& neighborId : neighbors) {
            auto neighbor = game.Get(neighborId);
            if (neighbor->alive) {
                if (auto road = std::dynamic_pointer_cast<Road>(neighbor)) {
                    if (!game.RoadConnections(id).contains(road->id)) {
                        game.AddConnection(id, road->id);
                        AddSegment(game, {dx, dy});
                        road->AddSegment(game, {-dx, -dy});
                    }
                } else if (auto source = std::dynamic_pointer_cast<Source>(neighbor)) {
                    game.AddConnection(id, source->id);
                    AddSegment(game, {dx, dy});
                } else if (auto target = std::dynamic_pointer_cast<Target>(neighbor)) {
                    game.AddConnection(id, target->id);
                    AddSegment(game, {dx, dy});
                }
            }
        }
    }
}

struct Car::S {
    S(Car& car) : car(car) {
    }

    virtual ~S() = default;

    virtual S* Update(Game& game) = 0;

    Car& car;
};

struct WaitingOf : public Car::S {
    WaitingOf(Car& car, S* next) : S(car), next(next) {
    }

    virtual S* Update(Game& game) override {
        --waiting;
        if (waiting > 0) {
            return this;
        }
        waiting = turns;
        return next;
    }

    S* next;
    int turns = 50;
    int waiting = turns;
};

struct Backward : public Car::S {
    Backward(Car& car, std::vector<GridXy> path) : S(car), path(std::move(path)) {
    }

    virtual S* Update(Game& game) override {
        if (pathInd < path.size()) {
            car.Move(path, pathInd, entering_, game);
            return this;
        }
        path.clear();
        pathInd = 0;
        car.speed = {};
        return nullptr;
    }

    std::vector<GridXy> path;
    std::size_t pathInd = 0;
    bool entering_ = false;
};

void Car::Move(std::vector<GridXy>& path, std::size_t& pathInd, bool& entering_, Game& game) {
    if (entering_) {
        auto gridXy = path[pathInd];
        if (GetCenter(shape) == game.FromGrid(gridXy) + SIDE_LENGTH / 2.f) {
            auto diff = path[pathInd + 1] - gridPosition;
            speed = {diff.x * Car::FAST_SPEED, diff.y * Car::FAST_SPEED};
            entering_ = false;
        }
    } else {
        if (speed == sf::Vector2f{0, 0}) {
            auto nextGridXy = path[1];
            auto diff = nextGridXy - gridPosition;
            speed = {diff.x * Car::FAST_SPEED, diff.y * Car::FAST_SPEED};
            return;
        }
        auto position = shape->getPosition();
        auto newPosition = position + speed;
        if (!game.IsIn(newPosition, gridPosition)) {
            gridPosition = game.ToGrid(newPosition);
            ++pathInd;
            if (pathInd + 1 < path.size()) {
                entering_ = true;
            }
        }
        if (pathInd + 1 == path.size() &&
            (Contains(static_cast<const Target&>(*game.Get(target)).shape, shape->getPosition()) ||
             Contains(static_cast<const Source&>(*game.Get(source)).shape, shape->getPosition())))
        {
            pathInd = path.size();
            speed = {};
        }
    }
    shape->move(speed);
}

struct Forward : public Car::S {
    Forward(Car& car, std::vector<GridXy> path) : S(car), path(std::move(path)) {
    }

    virtual S* Update(Game& game) override {
        if (pathInd < path.size()) {
            car.Move(path, pathInd, entering_, game);
            return this;
        }
        --std::static_pointer_cast<Target>(game.Get(car.target))->people;
        --std::static_pointer_cast<Target>(game.Get(car.target))->takenPeople;
        std::ranges::reverse(path);
        return new WaitingOf(car, new Backward(car, std::move(path)));
    }

    std::vector<GridXy> path;
    std::size_t pathInd = 0;
    bool entering_ = false;
};

void Car::Update(Game& game) {
    if (state_) {
        auto old = state_;
        state_ = state_->Update(game);
        if (old != state_) {
            delete old;
        }
    }
}

void Car::Go(auto carPath, Id targetId) {
    target = targetId;
    state_ = new Forward(*this, std::move(carPath));
}

#pragma clang diagnostic pop
