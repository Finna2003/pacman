﻿
#include <SFML/Graphics.hpp>
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <queue>
#include <functional> 
#include <utility>

using namespace std;

const int GRID_SIZE = 30;
const int GRID_WIDTH = 20;
const int GRID_HEIGHT = 16;

const int WIDTH = GRID_SIZE * GRID_WIDTH;
const int HEIGHT = GRID_SIZE * GRID_HEIGHT;

int pacmanX = 1;
int pacmanY = 1;

int ghost1X = 9;
int ghost1Y = 7;

int ghost2X = 14;
int ghost2Y = 8;

int ghost3X = 5;
int ghost3Y = 14;

int newX;
int newY;

int counter = 0;

std::vector<std::vector<int>> labyrinth = {
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1},
    {1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0,1},
    {1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,1},
    {1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,1},
    {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,1},
    {1, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0,1},
    {1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,1},
    {1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0,1},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1},
    {1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0,1},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1},
    {1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0,1},
    {1, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0,1},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1},
    {1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0,1},
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,1},
    {1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0,1},
    {1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0,0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1}
};


bool isGameOver = false;
void moveGhost(int pacmanX, int pacmanY, int& ghostX, int& ghostY) {
    if (pacmanX == ghostX && pacmanY == ghostY) {
        isGameOver = true;
    }

    int dx = pacmanX - ghostX;
    int dy = pacmanY - ghostY;

    int firstX = ghostX;
    int firstY = ghostY;

    int newX = ghostX;
    int newY = ghostY;

    if (std::abs(dx) >=std::abs(dy)) {
        if (dx > 0) {
            newX += 1;

            if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX==ghostX && newY==ghostY) ){
                if (dy > 0) {
                    newY += 1;
                    newX -= 1;
                    if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX == ghostX && newY == ghostY)) {
                        newY -= 2;

                    }
                }
                if (dy < 0) {
                    newY -= 1;

                    if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX == ghostX && newY == ghostY)) {
                        newY += 2;
                        newX -= 1;
                    }
                }
                else { newX -= 2;            
                }
            }
        }
        else if (dx < 0) {
            newX -= 1;
            if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX == ghostX && newY == ghostY)) {
                if (dy > 0) {
                    newY += 1;
                    newX += 1;
                }
                if (dy < 0) {
                    newY -= 1;
                    newX += 1;
                }
                else { newX += 2;             
                }
            }
        }
    }
    else if (std::abs(dx) < std::abs(dy)) {
        if (dy > 0) {
            newY += 1;

            if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX == ghostX && newY == ghostY)) {
                if (dx > 0) {
                    newX += 1;
                    newY -= 1;

                    if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX == ghostX && newY == ghostY)) {
                        newX -= 2;

                    }

                }
                if (dx < 0) {
                    newX -= 1;
                    newY -= 1;
                    if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX == ghostX && newY == ghostY)) {
                        newX += 2;

                    }
                }
                else { newY -= 2;             
                }
            }
        }
        else if (dy < 0) {
            newY -= 1;

            if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX == ghostX && newY == ghostY)) {
                if (dx > 0) {
                    newY += 1;
                    newX += 1;
                    if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX == ghostX && newY == ghostY)) {
                        newX -= 2;
                    }
                }
                if (dx < 0) {
                    newY += 1;
                    newX -= 1;

                    if (newX <= 0 || newX > GRID_WIDTH || newY <= 0 || newY > GRID_HEIGHT || labyrinth[newX][newY] != 0 || (newX == ghostX && newY == ghostY)) {
                        newX += 2;

                    }
                }
                else { newY += 2;            
                }
            }
        }
    }
    if (newX >= 0 && newX < GRID_WIDTH && newY >= 0 && newY < GRID_HEIGHT &&
        labyrinth[newX][newY] == 0) {
        ghostX = newX;
        ghostY = newY;
        std::cout << "Ghost moved to X: " << ghostX << ", Y: " << ghostY << std::endl;
    }

}
pair<int, int> start = make_pair(ghost3X, ghost3Y);

class Node {
public:
    pair<int, int> coordinates;
    Node* parent;
    double g, h;

    Node(pair<int, int> coord, Node* p, double gVal, double hVal)
        : coordinates(coord), parent(p), g(gVal), h(hVal) {}

    double getF() const { return g + h; }
};

bool isValid(int row, int col) {
    return (row >= 0) && (row < GRID_WIDTH) &&
        (col >= 0) && (col < GRID_HEIGHT);
}

bool isUnblocked(const vector<vector<int>>& maze, int row, int col) {
    return maze[row][col] == 0;
}

bool isDestination(int row, int col, pair<int, int> dest) {
    return (row == dest.first && col == dest.second);
}

double calculateHValue(int row, int col, pair<int, int> dest) {
    return sqrt((row - dest.first) * (row - dest.first) + (col - dest.second) * (col - dest.second));
}

void tracePath(Node* end) {
    if (end == nullptr)
        return;
    tracePath(end->parent);
    cout << "(" << end->coordinates.first << ", " << end->coordinates.second << ") " << endl;
}

void aStarSearch(const vector<vector<int>>& maze, pair<int, int> start, pair<int, int> dest) {
    if (!isValid(start.first, start.second) || !isValid(dest.first, dest.second)) {
        cout << "Invalid start or end\n";
        return;
    }

    if (!isUnblocked(maze, start.first, start.second) || !isUnblocked(maze, dest.first, dest.second)) {
        cout << "Start or end blocked\n";
        return;
    }

    if (isDestination(start.first, start.second, dest)) {
        cout << "We are already at our destination\n";
        return;
    }

    bool closedList[GRID_WIDTH][GRID_HEIGHT];
    memset(closedList, false, sizeof(closedList));

    priority_queue<Node*, vector<Node*>, function<bool(Node*, Node*)>> openList([](Node* a, Node* b) { return a->getF() > b->getF(); });

    openList.push(new Node(start, nullptr, 0.0, 0.0));

    vector<pair<int, int>> path;
    while (!openList.empty()) {
        Node* current = openList.top();
        openList.pop();


        int row = current->coordinates.first;
        int col = current->coordinates.second;
        closedList[row][col] = true;

        path.push_back(make_pair(row, col));
        counter += 1;

        if (isDestination(row, col, dest)) {
            cout << "Found a way: ";
            tracePath(current);
            return;
        }
        if (!path.empty()) {
            ghost3X = path[0].first;
            ghost3Y = path[0].second;
        }
        for (int i = -1; i <= 1; ++i) {
            for (int j = -1; j <= 1; ++j) {
                double gNew = current->g + 1.0;
                double hNew = calculateHValue(row + i, col + j, dest);
                double fNew = gNew + hNew;

                if (isValid(row + i, col + j) && isUnblocked(maze, row + i, col + j) && !closedList[row + i][col + j]) {
                    openList.push(new Node(make_pair(row + i, col + j), current, gNew, hNew));
                }
            }
        }
    }

    cout << "Path not found\n";
}


void movePacman(int dx, int dy) {
    if (isGameOver) {
        return;
    }

    int newX = pacmanX + dx;
    int newY = pacmanY + dy;

    if (newX >= 0 && newX < GRID_WIDTH && newY >= 0 && newY < GRID_HEIGHT &&
        labyrinth[newX][newY] == 0) {
        pacmanX = newX;
        pacmanY = newY;
    }


    if (pacmanX == ghost1X && pacmanY == ghost1Y ||
        pacmanX == ghost2X && pacmanY == ghost2Y ||
        pacmanX == ghost3X && pacmanY == ghost3Y) {
        isGameOver = true;

    }
    if (pacmanX == 18 && pacmanY == 15) {
        std::cout << "Win! " ;

    }

}


void drawGame(sf::RenderWindow& window) {
    window.clear();


    for (int x = 0; x < GRID_WIDTH; x++) {
        for (int y = 0; y < GRID_HEIGHT; y++) {
            if (labyrinth[x][y] == 1) {
                sf::RectangleShape wall(sf::Vector2f(GRID_SIZE, GRID_SIZE));
                wall.setPosition(x * GRID_SIZE, y * GRID_SIZE);
                wall.setFillColor(sf::Color::Yellow);
                window.draw(wall);
            }
        }
    }


    sf::RectangleShape pacman(sf::Vector2f(GRID_SIZE, GRID_SIZE));
    pacman.setPosition(pacmanX * GRID_SIZE, pacmanY * GRID_SIZE);
    pacman.setFillColor(sf::Color::Blue);
    window.draw(pacman);

    sf::RectangleShape ghost1(sf::Vector2f(GRID_SIZE, GRID_SIZE));
    ghost1.setPosition(ghost1X * GRID_SIZE, ghost1Y * GRID_SIZE);
    ghost1.setFillColor(sf::Color::Red);
    window.draw(ghost1);

    sf::RectangleShape ghost2(sf::Vector2f(GRID_SIZE, GRID_SIZE));
    ghost2.setPosition(ghost2X * GRID_SIZE, ghost2Y * GRID_SIZE);
    ghost2.setFillColor(sf::Color::Red);
    window.draw(ghost2);

    sf::RectangleShape ghost3(sf::Vector2f(GRID_SIZE, GRID_SIZE));
    ghost3.setPosition(ghost3X * GRID_SIZE, ghost3Y * GRID_SIZE);
    ghost3.setFillColor(sf::Color::Green);
    window.draw(ghost3);


    if (isGameOver) {
        return;
    }


    window.display();
}

int main() {
    sf::RenderWindow window(sf::VideoMode(WIDTH, HEIGHT), "Pacman Game");

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
            if (event.type == sf::Event::KeyPressed) {
                int first_x = pacmanX;
                int  first_y = pacmanY;
                if (event.key.code == sf::Keyboard::Up) {
                    movePacman(0, -1);
                }
                if (event.key.code == sf::Keyboard::Down) {
                    movePacman(0, 1);
                }
                if (event.key.code == sf::Keyboard::Left) {
                    movePacman(-1, 0);
                }
                if (event.key.code == sf::Keyboard::Right) {
                    movePacman(1, 0);
                }
                if (first_x != pacmanX || first_y != pacmanY) {
                    moveGhost(pacmanX, pacmanY, ghost1X, ghost1Y);
                    moveGhost(pacmanX, pacmanY, ghost2X, ghost2Y);
                    aStarSearch(labyrinth, start, make_pair(pacmanX, pacmanY));
                }
                if ((pacmanX == ghost1X && pacmanY == ghost1Y) ||(pacmanX == ghost3X && pacmanY == ghost3Y) ||(pacmanX == ghost2X && pacmanY == ghost2Y)){
                        std::cout << "Game over!" << std::endl;
                        return 0;
                    }

                
            }
        }

        

        drawGame(window);
    }
    return 0;
}


/*struct Point {
    double x, y;
};
int distSq(Point p1, Point p2) {
    int dx = p1.x - p2.x;
    int dy = p1.y - p2.y;
    return dx * dx + dy * dy;
}
Point pivot; 

int orientation(Point p, Point q, Point r) {
    int val = (q.y - p.y) * (r.x - q.x) - (q.x - p.x) * (r.y - q.y);
    if (val == 0) return 0;
    return (val > 0) ? 1 : -1;
}

bool compare(Point p, Point q) {
    int orientation_val = orientation(pivot, p, q);
    if (orientation_val == 0) {
        return (distSq(pivot, p) < distSq(pivot, q));
    }
    return (orientation_val == -1);
}

std::vector<Point> computeConvexHull(std::vector<Point>& points) {
    int min_x = points[0].x;
    int min_idx = 0;
    int n = points.size();
    for (int i = 1; i < n; i++) {
        int curr_x = points[i].x;
        if (curr_x < min_x || (curr_x == min_x && points[i].y < points[min_idx].y)) {
            min_x = curr_x;
            min_idx = i;
        }
    }

    std::swap(points[0], points[min_idx]);

    pivot = points[0];

    std::sort(points.begin() + 1, points.end(), compare);

    std::vector<Point> convexHull;
    convexHull.push_back(points[0]);
    convexHull.push_back(points[1]);
    convexHull.push_back(points[2]);

    for (int i = 3; i < n; i++) {
        while (convexHull.size() >= 2 && orientation(convexHull[convexHull.size() - 2], convexHull[convexHull.size() - 1], points[i]) != -1) {
            convexHull.pop_back();
        }
        convexHull.push_back(points[i]);
    }

    return convexHull;
}


double distance(const Point& p1, const Point& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

double getCircleRadius(const Point& center, const std::vector<Point>& points) {
    double radius = 0.0;
    for (const auto& point : points) {
        double dist = distance(center, point);
        radius = std::max(radius, dist);
    }
    return radius;
}

Point getCircleCenter(const Point& p1, const Point& p2) {
    double x = (p1.x + p2.x) / 2.0;
    double y = (p1.y + p2.y) / 2.0;
    return { x, y };
}

std::pair<Point, double> findMinimumEnclosingCircle(std::vector<Point>& points) {
    std::shuffle(points.begin(), points.end(), std::default_random_engine());
    Point center = points[0];
    double radius = 0.0;

    for (int i = 1; i < points.size(); ++i) {
        if (distance(center, points[i]) > radius) {
            center = points[i];
            radius = 0.0;

            for (int j = 0; j < i; ++j) {
                if (distance(center, points[j]) > radius) {
                    center = getCircleCenter(points[i], points[j]);
                    radius = getCircleRadius(center, std::vector<Point>(points.begin(), points.begin() + i + 1));
                }
            }
        }
    }

    return { center, radius };
}


int main() {
    int n;
    std::cout << "Enter count of points: ";
    std::cin >> n;

    std::vector<Point> points(n);
    std::cout << "Enter coordinates:\n";
    for (int i = 0; i < n; ++i) {
        std::cout << "Point  " << (i + 1) << " (x y): ";
        std::cin >> points[i].x >> points[i].y;
    }

    std::vector<Point> convexHull = computeConvexHull(points);

    std::cout << "Convex Hull Points:" << std::endl;
    for (const auto& point : convexHull) {
        std::cout << "(" << point.x << ", " << point.y << ")" << std::endl;
    }

    std::pair<Point, double> minimumEnclosingCircle = findMinimumEnclosingCircle(convexHull);

    Point center = minimumEnclosingCircle.first;
    double radius = minimumEnclosingCircle.second;

    std::cout << "Center: (" << center.x << ", " << center.y << ")" << std::endl;
    std::cout << "Radius: " << radius << std::endl;


    sf::RenderWindow window(sf::VideoMode(1200, 1000), "Circle");

    sf::CircleShape point(3); 
    point.setFillColor(sf::Color::Red); 

    sf::Vector2f centerPoint(center.x * 100, center.y * 100); // Центр кола
    float circleRadius = radius * 100; // Радіус кола


    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(sf::Color::White);

        sf::CircleShape point(3);
        point.setFillColor(sf::Color::Red);
        for (const auto& p : points) {
            point.setPosition(p.x * 100.f, p.y * 100.f);
            window.draw(point);
        }


        sf::CircleShape circle;
        circle.setRadius(circleRadius); 
        circle.setFillColor(sf::Color(0, 0, 255, 128));
        circle.setOutlineColor(sf::Color(0, 0, 255, 128)); 
        circle.setOutlineThickness(7); 

        // Задаємо позицію центра кола
        float centerX = centerPoint.x;
        float centerY = centerPoint.y;
        circle.setPosition(centerX - circle.getRadius(), centerY - circle.getRadius());
        window.draw(circle);

        window.display();
    }
    return 0;
}
*/


