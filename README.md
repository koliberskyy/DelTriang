# DelTriang
Триангуляция Делоне
Есть n точек на плоскости, проверить, что среди них есть выпуклый четырехугольник

Точек n <= 2 * 10^5
Координаты точек -1e9 <= x_i, y_i <= 1e9

Тип координат в задаче: числа с плавающей точкой. Решение будет тестироваться автоматически. Ожидается асимптитическая сложность (по времени и по памяти). 

Просьба присылать решение в одном файле main.cpp
решение должно считывать данные через std::cin
выводить через std::cout

Вывод в формате Yes или No

По времени иполнения, чтобы решение на 200000 точках работало не более чем за 2 секунды

подсказка: Create a Delauny triangulation. Find two adjacent triangles inside the convex hull (it means at least one of them is not located over convex hull) and report! To find more about Delauny triangulation read here.

As if there is any concave angle between four adjacent points, it will be flipped we can say every four adjacent points inside the convex hull create a convex quad.


@mfnbby->tg  koliberskyy@yandex.ru 
thx https://github.com/jbegaint/delaunay-cpp for base of this app
