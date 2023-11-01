#include "SStree.h"
#include <algorithm>
#include <utility>

bool SsNode::test(bool isRoot) const
{
    size_t count = 0;
    if (this->isLeaf())
    {
        const SsLeaf *leaf = dynamic_cast<const SsLeaf *>(this);
        count = leaf->points.size();

        // Verificar si los puntos están dentro del radio del nodo
        for (const Point &point : leaf->points)
        {
            if (distance(this->centroid, point) > this->radius)
            {
                std::cout << "Point outside node radius detected." << std::endl;
                return false;
            }
        }
    }
    else
    {
        const SsInnerNode *inner = dynamic_cast<const SsInnerNode *>(this);
        count = inner->children.size();

        // Verificar si los centroides de los hijos están dentro del radio del nodo padre
        for (const SsNode *child : inner->children)
        {
            if (distance(this->centroid, child->centroid) > this->radius)
            {
                std::cout << "Child centroid outside parent radius detected." << std::endl;
                return false;
            }
            // Verificar recursivamente cada hijo
            if (!child->test())
            {
                return false;
            }
        }
    }

    // Comprobar la validez de la cantidad de hijos/puntos
    if (!isRoot && (count < Settings::m || count > Settings::M))
    {
        std::cout << "Invalid number of children/points detected." << std::endl;
        return false;
    }

    // Comprobar punteros de parentezco, salvo para el nodo raíz
    if (!isRoot && !parent)
    {
        std::cout << "Node without parent detected." << std::endl;
        return false;
    }

    return true;
}

void SsTree::test() const
{
    bool result = root->test();

    if (root->parent)
    {
        std::cout << "Root node parent pointer is not null!" << std::endl;
        result = false;
    }

    if (result)
    {
        std::cout << "SS-Tree is valid!" << std::endl;
    }
    else
    {
        std::cout << "SS-Tree has issues!" << std::endl;
    }
}

void SsNode::print(size_t indent) const
{
    for (size_t i = 0; i < indent; ++i)
    {
        std::cout << "  ";
    }

    // Imprime información del nodo.
    std::cout << "Centroid: " << centroid << ", Radius: " << radius;
    if (isLeaf())
    {
        const SsLeaf *leaf = dynamic_cast<const SsLeaf *>(this);
        std::cout << ", Points: [ ";
        for (const Point &p : leaf->points)
        {
            std::cout << p << " ";
        }
        std::cout << "]";
    }
    else
    {
        std::cout << std::endl;
        const SsInnerNode *inner = dynamic_cast<const SsInnerNode *>(this);
        for (const SsNode *child : inner->children)
        {
            child->print(indent + 1);
        }
    }
    std::cout << std::endl;
}
void SsTree::print() const
{
    if (root)
    {
        root->print();
    }
    else
    {
        std::cout << "Empty tree." << std::endl;
    }
}

void SsLeaf::saveToStream(std::ostream &out) const
{
    // Guardar centroid
    auto D = centroid.dim();
    centroid.saveToFile(out, D);

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char *>(&radius_), sizeof(radius_));

    // Guardar el numero de puntos
    size_t numPoints = points.size();
    out.write(reinterpret_cast<const char *>(&numPoints), sizeof(numPoints));

    // Guardar los puntos
    for (const auto &point : points)
    {
        point.saveToFile(out, D);
    }

    // Guardar las rutas (paths)
    size_t numPaths = paths.size();
    out.write(reinterpret_cast<const char *>(&numPaths), sizeof(numPaths));
    for (const auto &p : paths)
    {
        size_t pathLength = p.size();
        out.write(reinterpret_cast<const char *>(&pathLength), sizeof(pathLength));
        out.write(p.c_str(), (long)pathLength);
    }
}

void SsInnerNode::saveToStream(std::ostream &out) const
{
    // Guardar centroid
    centroid.saveToFile(out, centroid.dim());

    // Guardar el radio
    float radius_ = radius.getValue();
    out.write(reinterpret_cast<const char *>(&radius_), sizeof(radius_));

    // Guardar si apunta a nodos hoja
    bool pointsToLeafs = children[0]->isLeaf();
    out.write(reinterpret_cast<const char *>(&pointsToLeafs), sizeof(pointsToLeafs));

    // Guardar la cantidad de hijos para saber cuántos nodos leer después
    size_t numChildren = children.size();
    out.write(reinterpret_cast<const char *>(&numChildren), sizeof(numChildren));

    // Guardar los hijos
    for (const auto &child : children)
    {
        child->saveToStream(out);
    }
}

void SsInnerNode::loadFromStream(std::istream &in)
{
    // Leer centroid
    auto D = centroid.dim();
    centroid.readFromFile(in, D);

    // leer el valor del radio
    float radius_ = 0;
    in.read(reinterpret_cast<char *>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // leer si apunta a hojas o nodos internos
    bool pointsToLeaf = false;
    in.read(reinterpret_cast<char *>(&pointsToLeaf), sizeof(pointsToLeaf));

    // leer cantidad de hijos
    size_t numChildren;
    in.read(reinterpret_cast<char *>(&numChildren), sizeof(numChildren));

    // leer hijos
    for (size_t i = 0; i < numChildren; ++i)
    {
        SsNode *child = pointsToLeaf ? static_cast<SsNode *>(new SsLeaf()) : static_cast<SsNode *>(new SsInnerNode());
        child->loadFromStream(in);
        children.push_back(child);
    }
}

void SsLeaf::loadFromStream(std::istream &in)
{
    // Leer centroid
    centroid.readFromFile(in, centroid.dim());

    // Leer radio
    float radius_ = 0;
    in.read(reinterpret_cast<char *>(&radius_), sizeof(radius_));
    this->radius = radius_;

    // Leer numero de puntos
    size_t numPoints;
    in.read(reinterpret_cast<char *>(&numPoints), sizeof(numPoints));

    // Leer puntos
    points.resize(numPoints);
    for (size_t i = 0; i < numPoints; ++i)
    {
        points[i].readFromFile(in, points[i].dim());
    }

    // Leer rutas (paths)
    size_t numPaths;
    in.read(reinterpret_cast<char *>(&numPaths), sizeof(numPaths));
    paths.resize(numPaths);
    for (size_t i = 0; i < numPaths; ++i)
    {
        size_t pathLength;
        in.read(reinterpret_cast<char *>(&pathLength), sizeof(pathLength));
        char *buffer = new char[pathLength + 1];
        in.read(buffer, (long)pathLength);
        buffer[pathLength] = '\0';
        paths[i] = std::string(buffer);
        delete[] buffer;
    }
}

void SsTree::saveToFile(const std::string &filename) const
{
    std::ofstream out(filename, std::ios::binary);
    if (!out)
    {
        throw std::runtime_error("Cannot open file for writing");
    }

    // Guardar las dimensiones de la estructura
    auto D = root->dim();
    out.write(reinterpret_cast<const char *>(&D), sizeof(D));

    // Guardar si el root es hija o nodo interno
    bool isLeaf = root->isLeaf();
    out.write(reinterpret_cast<const char *>(&isLeaf), sizeof(isLeaf));

    // Guardar el resto de la estructura
    root->saveToStream(out);
    out.close();
}

void SsTree::loadFromFile(const std::string &filename, size_t D)
{
    std::ifstream in(filename, std::ios::binary);
    if (!in)
    {
        throw std::runtime_error("Cannot open file for reading");
    }
    if (root)
    {
        delete root;
        root = nullptr;
    }

    // Aquí se asume que el primer valor determina las dimensiones
    in.read(reinterpret_cast<char *>(&D), sizeof(D));

    // El segundo valor determina si el root es hoja
    bool isLeaf;
    in.read(reinterpret_cast<char *>(&isLeaf), sizeof(isLeaf));
    if (isLeaf)
    {
        root = new SsLeaf();
    }
    else
    {
        root = new SsInnerNode();
    }
    root->loadFromStream(in);
    in.close();
}

void SsTree::insert(const Point &point)
{
    if (!root)
    {
        root = new SsLeaf();
        root->centroid = point;
        root->parent = nullptr;
        root->radius = 0;
        root->updateBoundingEnvelope();
        return;
    }

    std::pair<SsNode *, SsNode *> nodes = root->insert(point);
    if (nodes.first != nullptr)
    {
        root = new SsInnerNode();
        dynamic_cast<SsInnerNode *>(root)->children.push_back(nodes.first);
        nodes.first->parent = root;
        dynamic_cast<SsInnerNode *>(root)->children.push_back(nodes.second);
        nodes.second->parent = root;
        root->updateBoundingEnvelope();
        root->parent = nullptr;
    }
}

void SsTree::insert(Point &point, const std::string &path)
{
    point.path = "//home/renatoseb/2023-2/eda/labs/ss-tree/code/ss-tree/" + path;
    if (!root)
    {
        root = new SsLeaf();
        dynamic_cast<SsLeaf *>(root)->points.push_back(point);
        root->updateBoundingEnvelope();
        root->parent = nullptr;
    }
    else
    {
        std::pair<SsNode *, SsNode *> newChilds = root->insert(point);
        if (newChilds.first != nullptr)
        {
            root = new SsInnerNode();
            dynamic_cast<SsInnerNode *>(root)->children.push_back(newChilds.first);
            newChilds.first->parent = root;
            dynamic_cast<SsInnerNode *>(root)->children.push_back(newChilds.second);
            newChilds.second->parent = root;
            dynamic_cast<SsInnerNode *>(root)->updateBoundingEnvelope();
            dynamic_cast<SsInnerNode *>(root)->parent = nullptr;
        }
    }
}

std::pair<SsNode *, SsNode *> SsInnerNode::insert(const Point &point)
{
    SsNode *child = findClosestChild(point);
    std::pair<SsNode *, SsNode *> nodes = child->insert(point);
    if (nodes.first == nullptr)
    {
        this->updateBoundingEnvelope();
        return std::make_pair(nullptr, nullptr);
    }
    else
    {
        children.erase(remove(children.begin(), children.end(), child), children.end());
        nodes.first->parent = this;
        nodes.second->parent = this;
        children.push_back(nodes.first);
        children.push_back(nodes.second);
        updateBoundingEnvelope();
        if (children.size() > Settings::M)
        {
            return split();
        }
        else
        {
            return std::make_pair(nullptr, nullptr);
        }
    }
}

std::pair<SsNode *, SsNode *> SsLeaf::insert(const Point &point)
{
    if (std::find(points.begin(), points.end(), point) != points.end())
    {
        return std::make_pair(nullptr, nullptr);
    }
    points.push_back(point);
    updateBoundingEnvelope();
    if (points.size() > Settings::M)
    {
        return split();
    }
    return std::make_pair(nullptr, nullptr);
}

// NOTE: IMPLEMENTED AS BOOK
SsNode *SsInnerNode::findClosestChild(const Point &target) const
{
    SsNode *closest = nullptr;
    NType minDistance = std::numeric_limits<NType>::max();
    for (SsNode *child : children)
    {
        NType distance = ::distance(child->centroid, target);
        if (distance < minDistance)
        {
            minDistance = distance;
            closest = child;
        }
    }
    return closest;
}

// NOTE: IMPLEMENTED AS BOOK
SsNode *SsTree::searchParentLeaf(SsNode *node, const Point &target)
{
    if (node->isLeaf())
        return node;
    SsInnerNode *inner = dynamic_cast<SsInnerNode *>(node);
    SsNode *child = inner->findClosestChild(target);
    return searchParentLeaf(child, target);
}

// NOTE: IMPLEMENTED AS BOOK
void SsInnerNode::updateBoundingEnvelope()
{
    std::vector<Point> points = getEntriesCentroids();
    centroid = Point(children[0]->centroid.dim());
    for (SsNode *child : children)
    {
        for (size_t i = 0; i < centroid.dim(); ++i)
        {
            centroid[i] += child->centroid[i];
        }
    }
    for (size_t i = 0; i < centroid.dim(); ++i)
    {
        centroid[i] /= children.size();
    }
    radius = 0;
    for (SsNode *child : children)
    {
        radius = NType::max(radius, ::distance(centroid, child->centroid) + child->radius);
    }
}

// NOTE: IDK IF THIS IS OK, BUT IS THE SAME OF THE BOOK
void SsLeaf::updateBoundingEnvelope()
{
    std::vector<Point> points = getEntriesCentroids();
    centroid = Point(points[0].dim());
    for (const Point &point : points)
    {
        for (size_t i = 0; i < centroid.dim(); ++i)
        {
            centroid[i] += point[i];
        }
    }
    for (size_t i = 0; i < centroid.dim(); ++i)
    {
        centroid[i] /= points.size();
    }
    radius = 0;
    for (const Point &point : points)
    {
        radius = NType::max(radius, ::distance(centroid, point));
    }
}

// NOTE: NOT IMPLEMENTED IN BOOK, I SUPPOSE IS OK
NType SsNode::varianceAlongDirection(const std::vector<Point> &centroids, size_t direction) const
{
    NType mean = 0;
    for (const auto &centroid : centroids)
    {
        mean += centroid[direction];
    }
    mean /= centroids.size();

    NType variance = 0;
    for (const auto &centroid : centroids)
    {
        variance += NType::pow(centroid[direction] - mean, 2);
    }
    variance /= centroids.size();

    return variance;
}

// NOTE: IMPLEMENTED AS BOOK
size_t SsNode::directionOfMaxVariance() const
{
    NType maxVariance = 0;
    int directionIndex = 0;
    for (size_t i = 0; i < centroid.dim(); ++i)
    {
        NType variance = varianceAlongDirection(getEntriesCentroids(), i);
        if (variance > maxVariance)
        {
            maxVariance = variance;
            directionIndex = i;
        }
    }
    return directionIndex;
}

// NOTE: UNIQUE DIFFERENCE FROM BOOK (SORT CHILDREN INSTEAD OF CENTROIDS)
void SsInnerNode::sortEntriesByCoordinate(size_t coordinateIndex)
{
    std::sort(children.begin(), children.end(), [coordinateIndex](const SsNode *a, const SsNode *b)
              { return a->centroid[coordinateIndex] < b->centroid[coordinateIndex]; });
}

// NOTE: UNIQUE DIFFERENCE FROM BOOK (SORT POINTS INSTEAD OF CENTROIDS)
void SsLeaf::sortEntriesByCoordinate(size_t coordinateIndex)
{
    std::sort(points.begin(), points.end(), [coordinateIndex](const Point &a, const Point &b)
              { return a[coordinateIndex] < b[coordinateIndex]; });
}

// NOTE: IMPLEMENTED AS BOOK
size_t SsNode::minVarianceSplit(size_t coordinateIndex)
{
    NType minVariance = NType::INF;
    size_t splitIndex = Settings::m;
    for (size_t i = Settings::m; i < getEntriesCentroids().size() - Settings::m; ++i)
    {
        std::vector<Point> centroids1;
        std::vector<Point> centroids2;
        for (size_t j = 0; j < i - 1; ++j)
        {
            centroids1.push_back(getEntriesCentroids()[j]);
        }
        for (size_t j = i; j < getEntriesCentroids().size(); ++j)
        {
            centroids2.push_back(getEntriesCentroids()[j]);
        }
        // NOTE: UNIQUE DIFFERENCE FROM BOOK
        NType variance1 = varianceAlongDirection(centroids1, coordinateIndex);
        NType variance2 = varianceAlongDirection(centroids2, coordinateIndex);
        if (variance1 + variance2 < minVariance)
        {
            minVariance = variance1 + variance2;
            splitIndex = i;
        }
    }
    return splitIndex;
}

// NOTE: IMPLEMENTED AS BOOK
size_t SsNode::findSplitIndex()
{
    size_t index = directionOfMaxVariance();
    sortEntriesByCoordinate(index);
    // NOTE: UNIQUE DIFFERENCE FROM BOOK
    return minVarianceSplit(index);
}

// NOTE: IMPLEMENTED AS BOOK
std::pair<SsNode *, SsNode *> SsInnerNode::split()
{
    /*
    size_t splitIndex = findSplitIndex();
    SsInnerNode *left = new SsInnerNode();
    SsInnerNode *right = new SsInnerNode();
    for (size_t i = 0; i < splitIndex; ++i)
    {
        left->children.push_back(children[i]);
        children[i]->parent = left;
    }
    for (size_t i = splitIndex; i < children.size(); ++i)
    {
        right->children.push_back(children[i]);
        children[i]->parent = right;
    }

    left->updateBoundingEnvelope();
    right->updateBoundingEnvelope();

    return std::make_pair(left, right);
    */
    size_t splitIndex = findSplitIndex();
    std::vector<SsNode *> l, r;
    for (int i = 0; i < splitIndex + 1; i++)
    {
        l.push_back(children[i]);
    }
    for (int i = splitIndex + 1; i < children.size(); i++)
    {
        r.push_back(children[i]);
    }
    SsInnerNode *lnode = new SsInnerNode();
    SsInnerNode *rnode = new SsInnerNode();
    lnode->children = l;
    rnode->children = r;

    for (auto child : l)
    {
        child->parent = lnode;
    }

    for (auto child : r)
    {
        child->parent = rnode;
    }
    lnode->updateBoundingEnvelope();
    rnode->updateBoundingEnvelope();

    return std::make_pair(lnode, rnode);
}

// NOTE: IMPLEMENTED AS BOOK
std::pair<SsNode *, SsNode *> SsLeaf::split()
{

    size_t splitIndex = findSplitIndex();
    SsLeaf *left = new SsLeaf();
    SsLeaf *right = new SsLeaf();
    for (size_t i = 0; i < splitIndex; ++i)
    {
        left->points.push_back(points[i]);
    }
    for (size_t i = splitIndex; i < points.size(); ++i)
    {
        right->points.push_back(points[i]);
    }

    left->updateBoundingEnvelope();
    right->updateBoundingEnvelope();

    return std::make_pair(left, right);
}

// NOTE: IMPLEMENTED AS BOOK
std::vector<Point> SsLeaf::getEntriesCentroids() const
{
    return this->points;
}

// NOTE: IMPLEMENTED AS BOOK
std::vector<Point> SsInnerNode::getEntriesCentroids() const
{
    std::vector<Point> centroids;
    for (const auto &child : children)
    {
        centroids.push_back(child->centroid);
    }
    return centroids;
}

void SsInnerNode::FNDFTrav(const Point &q, size_t, std::priority_queue<Pair, std::vector<Pair>, Comparator> &L, NType &Dk) const
{
    for (const auto &child : children)
    {
        if (child->intersectsPoint(q))
        {
            child->FNDFTrav(q, 0, L, Dk);
        }
    }
}

void SsLeaf::FNDFTrav(const Point &q, size_t k, std::priority_queue<Pair, std::vector<Pair>, Comparator> &L, NType &Dk) const
{
    for (const auto &point : points)
    {
        NType distance = ::distance(point, q);
        if (distance < Dk)
        {
            Pair pair(point, distance);
            L.push(pair);
            if (L.size() > k)
            {
                L.pop();
                Dk = L.top().distance;
            }
        }
    }
}

std::vector<Point> SsTree::kNNQuery(const Point &center, size_t k) const
{
    std::priority_queue<Pair, std::vector<Pair>, Comparator> L;
    NType Dk = std::numeric_limits<NType>::max();
    root->FNDFTrav(center, k, L, Dk);
    std::vector<Point> result;
    while (!L.empty())
    {
        result.push_back(L.top().point);
        L.pop();
    }
    return result;
}