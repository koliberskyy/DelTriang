#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <cmath>

constexpr double eps = 1e-4;

template <typename T>
struct Point {
  T x, y;

  Point() : x{0}, y{0} {}
  Point(T _x, T _y) : x{_x}, y{_y} {}

  template <typename U>
  Point(U _x, U _y) : x{static_cast<T>(_x)}, y{static_cast<T>(_y)}
  {
  }

  friend std::ostream& operator<<(std::ostream& os, const Point<T>& p)
  {
    os << "x=" << p.x << "  y=" << p.y;
    return os;
  }

  bool operator==(const Point<T>& other) const
  {
    return (other.x == x && other.y == y);
  }

  bool operator!=(const Point<T>& other) const { return !operator==(other); }

};


template <typename T>
struct Edge {
  using Node = Point<T>;
  Node p0, p1;

  Edge(Node const& _p0, Node const& _p1) : p0{_p0}, p1{_p1} {}
  Edge() : p0{0, 0}, p1{0, 0} {}

  friend std::ostream& operator<<(std::ostream& os, const Edge& e)
  {
    os << "p0: [" << e.p0 << " ] p1: [" << e.p1 << "]";
    return os;
  }

  bool operator==(const Edge& other) const
  {
    return ((other.p0 == p0 && other.p1 == p1) ||
            (other.p0 == p1 && other.p1 == p0));
  }
};

template <typename T>
struct Circle {
  T x, y, radius;
  Circle() = default;
};

template <typename T>
struct Triangle {
  using Node = Point<T>;
  Node p0, p1, p2;
  Edge<T> e0, e1, e2;
  Circle<T> circle;

  Triangle(const Node& _p0, const Node& _p1, const Node& _p2)
      : p0{_p0},
        p1{_p1},
        p2{_p2},
        e0{_p0, _p1},
        e1{_p1, _p2},
        e2{_p0, _p2},
        circle{}
  {
    const auto ax = p1.x - p0.x;
    const auto ay = p1.y - p0.y;
    const auto bx = p2.x - p0.x;
    const auto by = p2.y - p0.y;

    const auto m = p1.x * p1.x - p0.x * p0.x + p1.y * p1.y - p0.y * p0.y;
    const auto u = p2.x * p2.x - p0.x * p0.x + p2.y * p2.y - p0.y * p0.y;
    const auto s = 1. / (2. * (ax * by - ay * bx));

    circle.x = ((p2.y - p0.y) * m + (p0.y - p1.y) * u) * s;
    circle.y = ((p0.x - p2.x) * m + (p1.x - p0.x) * u) * s;

    const auto dx = p0.x - circle.x;
    const auto dy = p0.y - circle.y;
    circle.radius = dx * dx + dy * dy;
  }
};

template <typename T>
struct Delaunay {
  std::vector<Triangle<T>> triangles;
  std::vector<Edge<T>> edges;
  bool isHasConvexSquare{false};
};

template <typename T>
Point<T> notIncludePoint(const Edge<T> edge, const Triangle<T> triangle){
    if(triangle.p0 != edge.p0 && triangle.p0 != edge.p1)
        return triangle.p0;
    else if(triangle.p1 != edge.p0 && triangle.p0 != edge.p1)
        return triangle.p1;
    else if(triangle.p2 != edge.p0 && triangle.p0 != edge.p1)
        return triangle.p2;
    return Point<T>();

}
template<typename T>
bool are_crossing(T x1, T x2, T  y1, T  y2, T  x3, T  x4, T  y3, T  y4)
{
    const T
 x = ((x1 * y2 - x2 * y1) * (x4 - x3) - (x3 * y4 - x4 * y3) * (x2 - x1)) / ((y1 - y2) * (x4 - x3) - (y3 - y4) * (x2 - x1)),
 y = ((y3 - y4) * x - (x3 * y4 - x4 * y3)) / (x4 - x3);

    if ((((x1 <= x) && (x2 >= x) && (x3 <= x) && (x4 >= x)) || ((y1 <= y) && (y2 >= y) && (y3 <= y) && (y4 >= y))))
        return true;
    else
        return false;
}

template<typename T>
bool are_crossing(Edge<T> lhs, Edge<T> rhs)
{
    return are_crossing(lhs.p0.x, lhs.p1.x, lhs.p0.y, lhs.p1.y, rhs.p0.x, rhs.p1.x, rhs.p0.y, rhs.p1.y);
}

template <typename T>
bool isСonvex(const Triangle<T> lhs, const Triangle<T> rhs){
    Edge<T> crossEdge;
    Edge<T> secondEdge;

    if(lhs.e0 == rhs.e0){
        crossEdge = lhs.e0;
    }
    else if(lhs.e1 == rhs.e1){
        crossEdge = lhs.e1;
    }
    else if(lhs.e2 == rhs.e2){
        crossEdge = lhs.e2;
    }
    else{
        return false;
    }

    secondEdge.p0 = notIncludePoint(crossEdge, lhs);
    secondEdge.p1 = notIncludePoint(crossEdge, rhs);



    return are_crossing(crossEdge, secondEdge);
}


template <
    typename T,
    typename = typename std::enable_if<std::is_floating_point<T>::value>::type>
Delaunay<T> triangulate(const std::vector<Point<T>>& points)
{
  using Node = Point<T>;
  if (points.size() < 3) {
    return Delaunay<T>{};
  }
  auto xmin = points[0].x;
  auto xmax = xmin;
  auto ymin = points[0].y;
  auto ymax = ymin;
  for (auto const& pt : points) {
    xmin = std::min(xmin, pt.x);
    xmax = std::max(xmax, pt.x);
    ymin = std::min(ymin, pt.y);
    ymax = std::max(ymax, pt.y);
  }

  const auto dx = xmax - xmin;
  const auto dy = ymax - ymin;
  const auto dmax = std::max(dx, dy);
  const auto midx = (xmin + xmax) / static_cast<T>(2.);
  const auto midy = (ymin + ymax) / static_cast<T>(2.);

  auto d = Delaunay<T>{};

  const auto p0 = Node{midx - 20 * dmax, midy - dmax};
  const auto p1 = Node{midx, midy + 20 * dmax};
  const auto p2 = Node{midx + 20 * dmax, midy - dmax};
  d.triangles.emplace_back(Triangle<T>{p0, p1, p2});

  for (auto const& pt : points) {
    std::vector<Edge<T>> edges;
    std::vector<Triangle<T>> tmps;
    for (auto const& tri : d.triangles) {
      const auto dist = (tri.circle.x - pt.x) * (tri.circle.x - pt.x) +
                        (tri.circle.y - pt.y) * (tri.circle.y - pt.y);
      if ((dist - tri.circle.radius) <= eps) {
        edges.push_back(tri.e0);
        edges.push_back(tri.e1);
        edges.push_back(tri.e2);
      }
      else {
        tmps.push_back(tri);
      }
    }

    std::vector<bool> remove(edges.size(), false);
    for (auto it1 = edges.begin(); it1 != edges.end(); ++it1) {
      for (auto it2 = edges.begin(); it2 != edges.end(); ++it2) {
        if (it1 == it2) {
          continue;
        }
        if (*it1 == *it2) {
          remove[std::distance(edges.begin(), it1)] = true;
          remove[std::distance(edges.begin(), it2)] = true;
        }
      }
    }

    edges.erase(
        std::remove_if(edges.begin(), edges.end(),
                       [&](auto const& e) { return remove[&e - &edges[0]]; }),
        edges.end());

    for (auto const& e : edges) {
      tmps.push_back({e.p0, e.p1, {pt.x, pt.y}});
    }

    //break triangulate
    if(tmps.size() > 1){
        auto existTriangles = tmps;
        existTriangles.erase(
                    std::remove_if(existTriangles.begin(), existTriangles.end(),
                                   [&](auto const& tri) {
                                     return ((tri.p0 == p0 || tri.p1 == p0 || tri.p2 == p0) ||
                                             (tri.p0 == p1 || tri.p1 == p1 || tri.p2 == p1) ||
                                             (tri.p0 == p2 || tri.p1 == p2 || tri.p2 == p2));
                                   }),
                    existTriangles.end());
        if(existTriangles.size() > 1){
            bool breakTrigger = false;

            for (auto const& lhs : existTriangles) {
                for (auto const& rhs : existTriangles) {
                    breakTrigger = isСonvex<T>(lhs, rhs);
                    if(breakTrigger){
                        d.isHasConvexSquare = true;
                        break;
                    }
                }
            }

            if(breakTrigger){
                d.triangles = tmps;
                break;
            }
        }
    }
    //break triangulate end

    d.triangles = tmps;
  }

  d.triangles.erase(
      std::remove_if(d.triangles.begin(), d.triangles.end(),
                     [&](auto const& tri) {
                       return ((tri.p0 == p0 || tri.p1 == p0 || tri.p2 == p0) ||
                               (tri.p0 == p1 || tri.p1 == p1 || tri.p2 == p1) ||
                               (tri.p0 == p2 || tri.p1 == p2 || tri.p2 == p2));
                     }),
      d.triangles.end());

  for (auto const& tri : d.triangles) {
    d.edges.push_back(tri.e0);
    d.edges.push_back(tri.e1);
    d.edges.push_back(tri.e2);
  }
  return d;
}


//формат ввода x_i, y_i, x_i, y_i ...
template<typename T>
std::vector<Point<T> > parceInput(const std::string& input){

    std::vector<Point<T>>  reply;
    auto end =  input.begin();
    auto begin = input.begin();
    while (end != input.end()) {
        end = std::find(begin, input.end(), ',');
        std::string tmpX(begin, end);

        if(end == input.end())
            break;
        begin = ++end;
        if( *begin == ' '){
            ++begin;
        }
        end = std::find(begin, input.end(), ',');
        std::string tmpY(begin, end);

        reply.emplace_back(Point<T>(atof(tmpX.c_str()), atof(tmpY.c_str())));

        if(end == input.end())
            break;
        begin = ++end;
        if( *begin == ' '){
            ++begin;
        }

    }

    return reply;
}
template<typename T>
std::string checkSquad(const Delaunay<T> &triangulation){
    if(triangulation.isHasConvexSquare){
        return "Yes";
    }
    return "No";
}

void test(const std::vector<Point<double> > &points){
    std::cout << "Input size:" << points.size() << '\n';
    auto start = std::chrono::system_clock::now();

    const auto triangulation = triangulate(points);
    std::cout << "Done: "<< checkSquad(triangulation) << '\n';

    auto timeCost = std::chrono::system_clock::now() - start;
    std::cout << "time cost:" << std::chrono::duration_cast<std::chrono::milliseconds>(timeCost).count() << " ms\n" ;
    std::cout << "mem cost:" << points.size() * sizeof (Point<double>) / 1000 << " kBytes\n--------------------------------- \n";
}

void randomTest(int testSize){
    std::vector<Point<double> > points;
    points.reserve(testSize);
    auto max = pow(10, 9);
    for(auto i = 0; i < testSize; i++){
        auto randomX = rand();
        auto randomY = rand();
        points.push_back(Point<double>(max / randomX, max/randomY));
    }
    test(points);
}

void oneLineTest(int testSize){

    std::cout << "one line test\n";
    std::vector<Point<double> > points;
    points.reserve(testSize);
    for(auto i = 0; i < testSize; i++){
        points.push_back(Point<double>(1.0, i));
    }
    test(points);
}

void oneLineOtherLastPointTest(int testSize){

    std::cout << "one line, other last point test\n";
    std::vector<Point<double> > points;
    points.reserve(testSize);
    for(auto i = 0; i < testSize; i++){
        points.push_back(Point<double>(1.0, i));
    }
    points.push_back(Point<double>(100.0, 10000));
    points.push_back(Point<double>(-500.1, 115000));

    test(points);
}

int main()
{
    //fake input
    test(std::vector<Point <double>>{
             Point<double>(1.0, 1.0),
             Point<double>(1.0, 1.0),
             Point<double>(1.0, 1.0),
             Point<double>(1.0, 1.0)
         });
    //fake input
    test(std::vector<Point <double>>{
             Point<double>(1.0, 10.0),
             Point<double>(2.0, 11.0),
             Point<double>(1.0, 1.0)
         });
    //fake input
     test(std::vector<Point <double>>{
              Point<double>(4.0, 1.0),
              Point<double>(3.0, 3.0),
              Point<double>(1.0, 3.0),
              Point<double>(3.0, 4.0)
          });

    //random inputs
    randomTest(1000);
    randomTest(10000);
    randomTest(100000);
    randomTest(200000);
    //one line tests
    oneLineTest(1000);
    oneLineOtherLastPointTest(1000);
    oneLineTest(10000);
    oneLineOtherLastPointTest(10000);
    oneLineTest(100000);
    oneLineTest(200000);

    std::string inputStr;
    getline(std::cin, inputStr);
    auto points = parceInput<double>(inputStr);
    const auto triangulation = triangulate(points);
    std::cout << checkSquad(triangulation);

    return 0;
}



