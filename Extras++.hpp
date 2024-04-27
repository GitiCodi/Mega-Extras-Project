#pragma once
// Includes

#include <iostream>
#include <cmath>
#include <limits>
#include <random>
#include <iomanip>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <list>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <stack>
#include <stdexcept>
#include <Windows.h>
#include <regex>
#include <cassert>
#include <memory>
#include <string>
#include <TlHelp32.h>
#include <psapi.h>
#include <math.h>
#include <numbers>
#include <functional>
#include <ctime>
#include <thread>
#include <vector>
#include <array>
#include <map>

// Structs

typedef struct Vector2
{
	int x, y;
};

typedef struct Vector3
{
	int x, y, z;
};

typedef struct Vector4
{
	int r, g, b, a;
};

// Program

class C_math
{
public:
	template <typename T>
	void add(T x, T y, T& resultVar) { resultVar = x + y; };

	template <typename T>
	void subtract(T x, T y, T& resultVar) { resultVar = x - y; };

	template <typename T>
	void multiply(T x, T y, T& resultVar) { resultVar = x * y; };

	template <typename T>
	void divide(T x, T y, T& resultVar) { resultVar = x / y; };
	
	template <typename T>
	void AreaRectangle(T width, T height, T& resultVar) { resultVar = width * height; };

	template <typename T>
	void AreaCuboid(T length, T& resultVar) { resultVar = length * length; };

	template <typename T>
	void AreaTriangle(T length, T height, T& resultVar) { resultVar = (length * height) / 2; };
};

using std::string;

class Memory
{
public:
    int GetProcessId(char* processName);
    int GetModuleBase(HANDLE processHandle, string& sModuleName);
    BOOL SetPrivilege(HANDLE hToken, LPCTSTR lpszPrivilege, BOOL bEnablePrivilege);
    BOOL GetDebugPrivileges(void);
    int ReadInt(HANDLE processHandle, int address);
    int GetPointerAddress(HANDLE processHandle, int startAddress, int offsets[], int offsetCount);
    int ReadPointerInt(HANDLE processHandle, int startAddress, int offsets[], int offsetCount);
    float ReadFloat(HANDLE processHandle, int address);
    float ReadPointerFloat(HANDLE processHandle, int startAddress, int offsets[], int offsetCount);
    char* ReadText(HANDLE processHandle, int address);
    char* ReadPointerText(HANDLE processHandle, int startAddress, int offsets[], int offsetCount);
};

int Memory::GetProcessId(char* processName) {
    SetLastError(0);
    PROCESSENTRY32 pe32;
    HANDLE hSnapshot = NULL;
    GetLastError();
    pe32.dwSize = sizeof(PROCESSENTRY32);
    hSnapshot = CreateToolhelp32Snapshot(TH32CS_SNAPPROCESS, 0);

    if (Process32First(hSnapshot, &pe32)) {
        do {
            if (strcmp((const char*)pe32.szExeFile, processName) == 0)
                break;
        } while (Process32Next(hSnapshot, &pe32));
    }

    if (hSnapshot != INVALID_HANDLE_VALUE)
        CloseHandle(hSnapshot);
    int err = GetLastError();
    //std::cout << err << std::endl;
    if (err != 0)
        return 0;
    return pe32.th32ProcessID;
}
int Memory::GetModuleBase(HANDLE processHandle, string& sModuleName)
{
    HMODULE* hModules = NULL;
    char szBuf[50];
    DWORD cModules;
    DWORD dwBase = -1;

    EnumProcessModules(processHandle, hModules, 0, &cModules);
    hModules = new HMODULE[cModules / sizeof(HMODULE)];

    if (EnumProcessModules(processHandle, hModules, cModules / sizeof(HMODULE), &cModules)) {
        for (size_t i = 0; i < cModules / sizeof(HMODULE); i++) {
            if (GetModuleBaseName(processHandle, hModules[i], (LPWSTR)szBuf, sizeof(szBuf))) {
                if (sModuleName.compare(szBuf) == 0) {
                    dwBase = (DWORD)hModules[i];
                    break;
                }
            }
        }
    }

    delete[] hModules;
    return dwBase;
}
BOOL Memory::SetPrivilege(HANDLE hToken, LPCTSTR lpszPrivilege, BOOL bEnablePrivilege)
{
    TOKEN_PRIVILEGES tp;
    LUID luid;

    if (!LookupPrivilegeValue(NULL, lpszPrivilege, &luid)) {
        //printf("LookupPrivilegeValue error: %u\n", GetLastError() );
        return FALSE;
    }

    tp.PrivilegeCount = 1;
    tp.Privileges[0].Luid = luid;
    if (bEnablePrivilege)
        tp.Privileges[0].Attributes = SE_PRIVILEGE_ENABLED;
    else
        tp.Privileges[0].Attributes = 0;

    if (!AdjustTokenPrivileges(hToken, FALSE, &tp, sizeof(TOKEN_PRIVILEGES), (PTOKEN_PRIVILEGES)NULL, (PDWORD)NULL)) {
        //printf("AdjustTokenPrivileges error: %u\n", GetLastError() );
        return FALSE;
    }

    if (GetLastError() == ERROR_NOT_ALL_ASSIGNED) {
        //printf("The token does not have the specified privilege. \n");
        return FALSE;
    }

    return TRUE;
}
BOOL Memory::GetDebugPrivileges(void) {
    HANDLE hToken = NULL;
    if (!OpenProcessToken(GetCurrentProcess(), TOKEN_ADJUST_PRIVILEGES, &hToken))
        return FALSE; //std::cout << "OpenProcessToken() failed, error\n>> " << GetLastError() << std::endl;
    //else std::cout << "OpenProcessToken() is OK, got the handle!" << std::endl;

    if (!SetPrivilege(hToken, SE_DEBUG_NAME, TRUE))
        return FALSE; //std::cout << "Failed to enable privilege, error:\n>> " << GetLastError() << std::endl;

    return TRUE;
}
int Memory::ReadInt(HANDLE processHandle, int address) {
    if (address == -1)
        return -1;
    int buffer = 0;
    SIZE_T NumberOfBytesToRead = sizeof(buffer); //this is equal to 4
    SIZE_T NumberOfBytesActuallyRead;
    BOOL success = ReadProcessMemory(processHandle, (LPCVOID)address, &buffer, NumberOfBytesToRead, &NumberOfBytesActuallyRead);
    if (!success || NumberOfBytesActuallyRead != NumberOfBytesToRead) {
        std::cout << "Memory Error!" << std::endl;
        return -1;
    }
    //if (err || NumberOfBytesActuallyRead != NumberOfBytesToRead) {
    //	DWORD lastError = GetLastError();
    //	if (lastError != 0)
    //        std::cout << lastError << std::endl;
    //    std::cout << "blub" << std::endl;
    //}
    return buffer;
}
int Memory::GetPointerAddress(HANDLE processHandle, int startAddress, int offsets[], int offsetCount) {
    if (startAddress == -1)
        return -1;
    int ptr = ReadInt(processHandle, startAddress);
    for (int i = 0; i < offsetCount - 1; i++) {
        ptr += offsets[i];
        ptr = ReadInt(processHandle, ptr);
    }
    ptr += offsets[offsetCount - 1];
    return ptr;
}
int Memory::ReadPointerInt(HANDLE processHandle, int startAddress, int offsets[], int offsetCount) {
    if (startAddress == -1)
        return -1;
    return ReadInt(processHandle, GetPointerAddress(processHandle, startAddress, offsets, offsetCount));
}
float Memory::ReadFloat(HANDLE processHandle, int address) {
    if (address == -1)
        return -1;
    float buffer = 0.0;
    SIZE_T NumberOfBytesToRead = sizeof(buffer); //this is equal to 4
    SIZE_T NumberOfBytesActuallyRead;
    BOOL success = ReadProcessMemory(processHandle, (LPCVOID)address, &buffer, NumberOfBytesToRead, &NumberOfBytesActuallyRead);
    if (!success || NumberOfBytesActuallyRead != NumberOfBytesToRead)
        return -1;
    return buffer;
}
float Memory::ReadPointerFloat(HANDLE processHandle, int startAddress, int offsets[], int offsetCount) {
    if (startAddress == -1)
        return -1;
    return ReadFloat(processHandle, GetPointerAddress(processHandle, startAddress, offsets, offsetCount));
}
char* Memory::ReadText(HANDLE processHandle, int address) {
    if (address == -1)
        return (char*) "-1";
    char buffer = !0;
    char* stringToRead = new char[128];
    SIZE_T NumberOfBytesToRead = sizeof(buffer);
    SIZE_T NumberOfBytesActuallyRead;
    int i = 0;
    while (buffer != 0) {
        BOOL success = ReadProcessMemory(processHandle, (LPCVOID)address, &buffer, NumberOfBytesToRead, &NumberOfBytesActuallyRead);
        if (!success || NumberOfBytesActuallyRead != NumberOfBytesToRead)
            return (char*) "-1";
        stringToRead[i] = buffer;
        i++;
        address++;
    }
    return stringToRead;
}
char* Memory::ReadPointerText(HANDLE processHandle, int startAddress, int offsets[], int offsetCount) {
    if (startAddress == -1)
        return (char*) "-1";
    return ReadText(processHandle, GetPointerAddress(processHandle, startAddress, offsets, offsetCount));
}

class Vector2D {
private:
    double x;
    double y;

public:
    // Constructors
    Vector2D() : x(0.0), y(0.0) {}
    Vector2D(double x_val, double y_val) : x(x_val), y(y_val) {}

    // Getters
    double getX() const { return x; }
    double getY() const { return y; }

    // Setters
    void setX(double x_val) { x = x_val; }
    void setY(double y_val) { y = y_val; }

    // Vector operations
    Vector2D operator+(const Vector2D& other) const {
        return Vector2D(x + other.x, y + other.y);
    }

    Vector2D operator-(const Vector2D& other) const {
        return Vector2D(x - other.x, y - other.y);
    }

    double dotProduct(const Vector2D& other) const {
        return x * other.x + y * other.y;
    }

    double crossProduct(const Vector2D& other) const {
        return x * other.y - y * other.x;
    }

    double magnitude() const {
        return sqrt(x * x + y * y);
    }

    // Normalize the vector
    Vector2D normalize() const {
        double mag = magnitude();
        return Vector2D(x / mag, y / mag);
    }
};

class Matrix {
private:
    std::vector<std::vector<double>> data;
    int rows;
    int cols;

public:
    // Constructors
    Matrix(int numRows, int numCols) : rows(numRows), cols(numCols) {
        data.resize(rows, std::vector<double>(cols, 0.0));
    }

    // Accessors
    int numRows() const { return rows; }
    int numCols() const { return cols; }

    // Element access
    double& operator()(int row, int col) {
        return data[row][col];
    }

    const double& operator()(int row, int col) const {
        return data[row][col];
    }

    // Matrix operations
    Matrix operator+(const Matrix& other) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] + other(i, j);
            }
        }
        return result;
    }

    Matrix operator-(const Matrix& other) const {
        Matrix result(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(i, j) = data[i][j] - other(i, j);
            }
        }
        return result;
    }

    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            throw std::invalid_argument("Matrix dimensions mismatch");
        }

        Matrix result(rows, other.cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                for (int k = 0; k < cols; ++k) {
                    result(i, j) += data[i][k] * other(k, j);
                }
            }
        }
        return result;
    }

    // Transpose the matrix
    Matrix transpose() const {
        Matrix result(cols, rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result(j, i) = data[i][j];
            }
        }
        return result;
    }
};

template <typename T>
class Node {
public:
    T data;
    Node* next;

    Node(const T& newData) : data(newData), next(nullptr) {}
};

template <typename T>
class LinkedList {
private:
    Node<T>* head;
    int size;

public:
    // Constructor
    LinkedList() : head(nullptr), size(0) {}

    // Destructor
    ~LinkedList() {
        clear();
    }

    // Insertion at the end
    void insert(const T& newData) {
        Node<T>* newNode = new Node<T>(newData);
        if (head == nullptr) {
            head = newNode;
        }
        else {
            Node<T>* current = head;
            while (current->next != nullptr) {
                current = current->next;
            }
            current->next = newNode;
        }
        size++;
    }

    // Deletion of the entire list
    void clear() {
        Node<T>* current = head;
        while (current != nullptr) {
            Node<T>* temp = current;
            current = current->next;
            delete temp;
        }
        head = nullptr;
        size = 0;
    }

    // Display the list
    void display() const {
        Node<T>* current = head;
        while (current != nullptr) {
            std::cout << current->data << " ";
            current = current->next;
        }
        std::cout << std::endl;
    }

    // Getter for size
    int getSize() const {
        return size;
    }
};

template <typename T>
class Stack {
private:
    static const int MAX_SIZE = 100;
    T data[MAX_SIZE];
    int topIndex;

public:
    // Constructor
    Stack() : topIndex(-1) {}

    // Push an element onto the stack
    void push(const T& element) {
        if (topIndex >= MAX_SIZE - 1) {
            throw std::overflow_error("Stack overflow");
        }
        data[++topIndex] = element;
    }

    // Pop an element from the stack
    T pop() {
        if (topIndex < 0) {
            throw std::underflow_error("Stack underflow");
        }
        return data[topIndex--];
    }

    // Peek at the top element without removing it
    T peek() const {
        if (topIndex < 0) {
            throw std::underflow_error("Stack is empty");
        }
        return data[topIndex];
    }

    // Check if the stack is empty
    bool isEmpty() const {
        return topIndex == -1;
    }

    // Get the number of elements in the stack
    int size() const {
        return topIndex + 1;
    }
};

template <typename T>
class Queue {
private:
    static const int MAX_SIZE = 100;
    T data[MAX_SIZE];
    int frontIndex;
    int rearIndex;

public:
    // Constructor
    Queue() : frontIndex(-1), rearIndex(-1) {}

    // Enqueue an element into the queue
    void enqueue(const T& element) {
        if ((rearIndex + 1) % MAX_SIZE == frontIndex) {
            throw std::overflow_error("Queue overflow");
        }
        if (frontIndex == -1) {
            frontIndex = 0;
        }
        rearIndex = (rearIndex + 1) % MAX_SIZE;
        data[rearIndex] = element;
    }

    // Dequeue an element from the queue
    T dequeue() {
        if (frontIndex == -1) {
            throw std::underflow_error("Queue underflow");
        }
        T element = data[frontIndex];
        if (frontIndex == rearIndex) {
            frontIndex = -1;
            rearIndex = -1;
        }
        else {
            frontIndex = (frontIndex + 1) % MAX_SIZE;
        }
        return element;
    }

    // Peek at the front element without removing it
    T peek() const {
        if (frontIndex == -1) {
            throw std::underflow_error("Queue is empty");
        }
        return data[frontIndex];
    }

    // Check if the queue is empty
    bool isEmpty() const {
        return frontIndex == -1;
    }

    // Get the number of elements in the queue
    int size() const {
        if (frontIndex == -1) {
            return 0;
        }
        return (rearIndex >= frontIndex) ? (rearIndex - frontIndex + 1) : (MAX_SIZE - frontIndex + rearIndex + 1);
    }
};

template <typename T>
class TreeNode {
public:
    T data;
    TreeNode* left;
    TreeNode* right;

    TreeNode(const T& value) : data(value), left(nullptr), right(nullptr) {}
};

template <typename T>
class BinarySearchTree {
private:
    TreeNode<T>* root;

    // Helper function to recursively insert a node into the BST
    void insertNode(TreeNode<T>*& node, const T& value) {
        if (node == nullptr) {
            node = new TreeNode<T>(value);
        }
        else {
            if (value < node->data) {
                insertNode(node->left, value);
            }
            else {
                insertNode(node->right, value);
            }
        }
    }

    // Helper function to recursively traverse the BST in inorder
    void inorderTraversal(TreeNode<T>* node, std::function<void(const T&)> visit) const {
        if (node != nullptr) {
            inorderTraversal(node->left, visit);
            visit(node->data);
            inorderTraversal(node->right, visit);
        }
    }

public:
    // Constructor
    BinarySearchTree() : root(nullptr) {}

    // Destructor
    ~BinarySearchTree() {
        clear(root);
    }

    // Insert an element into the BST
    void insert(const T& value) {
        insertNode(root, value);
    }

    // Inorder traversal of the BST
    void inorder(std::function<void(const T&)> visit) const {
        inorderTraversal(root, visit);
    }

    // Clear the BST
    void clear() {
        clear(root);
        root = nullptr;
    }

private:
    // Helper function to recursively clear the BST
    void clear(TreeNode<T>* node) {
        if (node != nullptr) {
            clear(node->left);
            clear(node->right);
            delete node;
        }
    }
};

template <typename T>
class Graph {
private:
    std::unordered_map<T, std::unordered_set<T>> adjacencyList;

public:
    // Add a vertex to the graph
    void addVertex(const T& vertex) {
        if (adjacencyList.find(vertex) == adjacencyList.end()) {
            adjacencyList[vertex] = std::unordered_set<T>();
        }
    }

    // Add an edge between two vertices
    void addEdge(const T& source, const T& destination) {
        addVertex(source);
        addVertex(destination);
        adjacencyList[source].insert(destination);
        // For undirected graph, uncomment the line below:
        // adjacencyList[destination].insert(source);
    }

    // Depth-first search traversal
    void dfs(const T& startVertex, std::function<void(const T&)> visit) const {
        std::unordered_set<T> visited;
        std::stack<T> stack;
        stack.push(startVertex);

        while (!stack.empty()) {
            T current = stack.top();
            stack.pop();

            if (visited.find(current) == visited.end()) {
                visit(current);
                visited.insert(current);

                for (const T& neighbor : adjacencyList[current]) {
                    if (visited.find(neighbor) == visited.end()) {
                        stack.push(neighbor);
                    }
                }
            }
        }
    }

    // Breadth-first search traversal
    void bfs(const T& startVertex, std::function<void(const T&)> visit) const {
        std::unordered_set<T> visited;
        std::queue<T> queue;
        queue.push(startVertex);

        while (!queue.empty()) {
            T current = queue.front();
            queue.pop();

            if (visited.find(current) == visited.end()) {
                visit(current);
                visited.insert(current);

                for (const T& neighbor : adjacencyList[current]) {
                    if (visited.find(neighbor) == visited.end()) {
                        queue.push(neighbor);
                    }
                }
            }
        }
    }
};

template<typename KeyType, typename ValueType>
class HashTable {
private:
    static const int TABLE_SIZE = 100;
    std::vector<std::list<std::pair<KeyType, ValueType>>> table;

    // Hash function
    size_t hash(const KeyType& key) const {
        return std::hash<KeyType>{}(key) % TABLE_SIZE;
    }

public:
    // Constructor
    HashTable() : table(TABLE_SIZE) {}

    // Insert key-value pair into the hash table
    void insert(const KeyType& key, const ValueType& value) {
        size_t index = hash(key);
        for (const auto& pair : table[index]) {
            if (pair.first == key) {
                throw std::runtime_error("Duplicate key");
            }
        }
        table[index].push_back(std::make_pair(key, value));
    }

    // Get the value associated with the given key
    ValueType get(const KeyType& key) const {
        size_t index = hash(key);
        for (const auto& pair : table[index]) {
            if (pair.first == key) {
                return pair.second;
            }
        }
        throw std::out_of_range("Key not found");
    }

    // Remove the key-value pair associated with the given key
    void remove(const KeyType& key) {
        size_t index = hash(key);
        auto& bucket = table[index];
        for (auto it = bucket.begin(); it != bucket.end(); ++it) {
            if (it->first == key) {
                bucket.erase(it);
                return;
            }
        }
        throw std::out_of_range("Key not found");
    }

    // Check if the hash table contains the given key
    bool contains(const KeyType& key) const {
        size_t index = hash(key);
        for (const auto& pair : table[index]) {
            if (pair.first == key) {
                return true;
            }
        }
        return false;
    }
};

class StringTokenizer {
private:
    std::vector<std::string> tokens;
    std::string delimiters;

public:
    // Constructor
    StringTokenizer(const std::string& str, const std::string& delims) : delimiters(delims) {
        tokenize(str);
    }

    // Tokenize the input string
    void tokenize(const std::string& str) {
        std::string token;
        std::istringstream iss(str);
        while (std::getline(iss, token, delimiters[0])) {
            tokens.push_back(token);
        }
    }

    // Get the tokens
    const std::vector<std::string>& getTokens() const {
        return tokens;
    }
};

namespace fs = std::filesystem;

class FileManager {
public:
    // Constructor
    FileManager() {}

    // Create a new file with the given content
    void createFile(const std::string& filename, const std::string& content) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << content;
            file.close();
            std::cout << "File \"" << filename << "\" created successfully." << std::endl;
        }
        else {
            std::cerr << "Error: Unable to create file \"" << filename << "\"" << std::endl;
        }
    }

    // Read the contents of a file into a string
    std::string readFile(const std::string& filename) const {
        std::ifstream file(filename);
        std::stringstream buffer;
        if (file.is_open()) {
            buffer << file.rdbuf();
            file.close();
        }
        else {
            std::cerr << "Error: Unable to read file \"" << filename << "\"" << std::endl;
        }
        return buffer.str();
    }

    // Copy a file to a new location
    void copyFile(const std::string& source, const std::string& destination) const {
        try {
            fs::copy_file(source, destination, fs::copy_options::overwrite_existing);
            std::cout << "File \"" << source << "\" copied to \"" << destination << "\" successfully." << std::endl;
        }
        catch (const fs::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }

    // Move a file to a new location
    void moveFile(const std::string& source, const std::string& destination) const {
        try {
            fs::rename(source, destination);
            std::cout << "File \"" << source << "\" moved to \"" << destination << "\" successfully." << std::endl;
        }
        catch (const fs::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }

    // Rename a file
    void renameFile(const std::string& oldName, const std::string& newName) const {
        try {
            fs::rename(oldName, newName);
            std::cout << "File \"" << oldName << "\" renamed to \"" << newName << "\" successfully." << std::endl;
        }
        catch (const fs::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }

    // Remove a file
    void removeFile(const std::string& filename) const {
        try {
            fs::remove(filename);
            std::cout << "File \"" << filename << "\" removed successfully." << std::endl;
        }
        catch (const fs::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
    }

    // Check if a file exists
    bool fileExists(const std::string& filename) const {
        return fs::exists(filename);
    }

    // List contents of a directory
    std::vector<std::string> listDirectory(const std::string& path) const {
        std::vector<std::string> contents;
        try {
            for (const auto& entry : fs::directory_iterator(path)) {
                contents.push_back(entry.path().filename().string());
            }
        }
        catch (const fs::filesystem_error& e) {
            std::cerr << "Error: " << e.what() << std::endl;
        }
        return contents;
    }
};


class Date {
private:
    int day;
    int month;
    int year;

public:
    // Constructor
    Date(int d, int m, int y) : day(d), month(m), year(y) {}

    // Getter methods
    int getDay() const { return day; }
    int getMonth() const { return month; }
    int getYear() const { return year; }

    // Setter methods
    void setDay(int d) { day = d; }
    void setMonth(int m) { month = m; }
    void setYear(int y) { year = y; }

    // Format date as string (DD/MM/YYYY)
    std::string toString() const {
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(2) << day << "/";
        ss << std::setfill('0') << std::setw(2) << month << "/";
        ss << year;
        return ss.str();
    }
};

class Time {
private:
    int hours;
    int minutes;
    int seconds;

public:
    // Constructor
    Time(int h, int m, int s) : hours(h), minutes(m), seconds(s) {}

    // Getter methods
    int getHours() const { return hours; }
    int getMinutes() const { return minutes; }
    int getSeconds() const { return seconds; }

    // Setter methods
    void setHours(int h) { hours = h; }
    void setMinutes(int m) { minutes = m; }
    void setSeconds(int s) { seconds = s; }

    // Format time as string (HH:MM:SS)
    std::string toString() const {
        std::stringstream ss;
        ss << std::setfill('0') << std::setw(2) << hours << ":";
        ss << std::setfill('0') << std::setw(2) << minutes << ":";
        ss << std::setfill('0') << std::setw(2) << seconds;
        return ss.str();
    }
};

class RandomNumberGenerator {
private:
    std::mt19937_64 engine; // Mersenne Twister 64-bit PRNG

public:
    // Constructor with seed
    RandomNumberGenerator(unsigned long long seed) : engine(seed) {}

    // Generate a random integer in the range [min, max]
    int generateInt(int min, int max) {
        std::uniform_int_distribution<int> distribution(min, max);
        return distribution(engine);
    }

    // Generate a random floating-point number in the range [min, max]
    double generateDouble(double min, double max) {
        std::uniform_real_distribution<double> distribution(min, max);
        return distribution(engine);
    }
};

class Color {
private:
    int red;
    int green;
    int blue;
    int alpha;

public:
    // Constructors
    Color() : red(0), green(0), blue(0), alpha(255) {}
    Color(int r, int g, int b, int a = 255) : red(r), green(g), blue(b), alpha(a) {}

    // Getter methods
    int getRed() const { return red; }
    int getGreen() const { return green; }
    int getBlue() const { return blue; }
    int getAlpha() const { return alpha; }

    // Setter methods
    void setRed(int r) { red = r; }
    void setGreen(int g) { green = g; }
    void setBlue(int b) { blue = b; }
    void setAlpha(int a) { alpha = a; }

    // Format color as string (R,G,B,A)
    std::string toString() const {
        std::stringstream ss;
        ss << "RGBA(" << red << "," << green << "," << blue << "," << alpha << ")";
        return ss.str();
    }
};

class MatrixMath {
public:
    // Add two matrices
    static std::vector<std::vector<double>> add(const std::vector<std::vector<double>>& matrix1,
        const std::vector<std::vector<double>>& matrix2) {
        if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
            throw std::invalid_argument("Matrix dimensions mismatch");
        }

        std::vector<std::vector<double>> result(matrix1.size(), std::vector<double>(matrix1[0].size(), 0.0));
        for (size_t i = 0; i < matrix1.size(); ++i) {
            for (size_t j = 0; j < matrix1[0].size(); ++j) {
                result[i][j] = matrix1[i][j] + matrix2[i][j];
            }
        }
        return result;
    }

    // Subtract one matrix from another
    static std::vector<std::vector<double>> subtract(const std::vector<std::vector<double>>& matrix1,
        const std::vector<std::vector<double>>& matrix2) {
        if (matrix1.size() != matrix2.size() || matrix1[0].size() != matrix2[0].size()) {
            throw std::invalid_argument("Matrix dimensions mismatch");
        }

        std::vector<std::vector<double>> result(matrix1.size(), std::vector<double>(matrix1[0].size(), 0.0));
        for (size_t i = 0; i < matrix1.size(); ++i) {
            for (size_t j = 0; j < matrix1[0].size(); ++j) {
                result[i][j] = matrix1[i][j] - matrix2[i][j];
            }
        }
        return result;
    }

    // Multiply two matrices
    static std::vector<std::vector<double>> multiply(const std::vector<std::vector<double>>& matrix1,
        const std::vector<std::vector<double>>& matrix2) {
        if (matrix1[0].size() != matrix2.size()) {
            throw std::invalid_argument("Matrix dimensions mismatch");
        }

        std::vector<std::vector<double>> result(matrix1.size(), std::vector<double>(matrix2[0].size(), 0.0));
        for (size_t i = 0; i < matrix1.size(); ++i) {
            for (size_t j = 0; j < matrix2[0].size(); ++j) {
                for (size_t k = 0; k < matrix2.size(); ++k) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
        return result;
    }
};

class CaesarCipher {
private:
    int shift;

public:
    // Constructor
    CaesarCipher(int s) : shift(s) {}

    // Encrypt plaintext using Caesar cipher
    std::string encrypt(const std::string& plaintext) const {
        std::string ciphertext;
        for (char c : plaintext) {
            if (std::isalpha(c)) {
                char base = std::isupper(c) ? 'A' : 'a';
                ciphertext += (c - base + shift) % 26 + base;
            }
            else {
                ciphertext += c;
            }
        }
        return ciphertext;
    }

    // Decrypt ciphertext using Caesar cipher
    std::string decrypt(const std::string& ciphertext) const {
        std::string plaintext;
        for (char c : ciphertext) {
            if (std::isalpha(c)) {
                char base = std::isupper(c) ? 'A' : 'a';
                plaintext += (c - base - shift + 26) % 26 + base;
            }
            else {
                plaintext += c;
            }
        }
        return plaintext;
    }
};

class Compression {
public:
    // Compress a string using run-length encoding
    static std::string compress(const std::string& input) {
        std::string compressed;
        size_t count = 1;
        for (size_t i = 0; i < input.size(); ++i) {
            if (i + 1 < input.size() && input[i] == input[i + 1]) {
                count++;
            }
            else {
                compressed += input[i];
                compressed += std::to_string(count);
                count = 1;
            }
        }
        return compressed;
    }

    // Decompress a string using run-length encoding
    static std::string decompress(const std::string& input) {
        std::string decompressed;
        for (size_t i = 0; i < input.size(); i += 2) {
            char ch = input[i];
            int count = std::stoi(input.substr(i + 1, 1));
            decompressed.append(count, ch);
        }
        return decompressed;
    }
};

class Geometry {
public:
    // Calculate the area of a circle with given radius
    static double circleArea(double radius) {
        return pi * radius * radius;
    }

    // Calculate the perimeter of a circle with given radius
    static double circlePerimeter(double radius) {
        return 2 * pi * radius;
    }

    // Calculate the area of a rectangle with given width and height
    static double rectangleArea(double width, double height) {
        return width * height;
    }

    // Calculate the perimeter of a rectangle with given width and height
    static double rectanglePerimeter(double width, double height) {
        return 2 * (width + height);
    }

    // Calculate the area of a triangle with given base and height
    static double triangleArea(double base, double height) {
        return 0.5 * base * height;
    }

    // Calculate the perimeter of a triangle with given side lengths
    static double trianglePerimeter(double side1, double side2, double side3) {
        return side1 + side2 + side3;
    }
public:
    static constexpr double pi = 3.14159265358979323846;
};

class Basic_Trigonometry {
public:
    // Constants
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double DEG_TO_RAD = PI / 180.0;
    static constexpr double RAD_TO_DEG = 180.0 / PI;

    // Basic trigonometric functions
    static double sin(double angle) { return std::sin(angle); }
    static double cos(double angle) { return std::cos(angle); }
    static double tan(double angle) { return std::tan(angle); }
    static double cot(double angle) { return 1.0 / std::tan(angle); }
    static double sec(double angle) { return 1.0 / std::cos(angle); }
    static double csc(double angle) { return 1.0 / std::sin(angle); }

    // Inverse trigonometric functions
    static double arcsin(double value) { return std::asin(value); }
    static double arccos(double value) { return std::acos(value); }
    static double arctan(double value) { return std::atan(value); }
    static double arccot(double value) { return std::atan(1.0 / value); }
    static double arcsec(double value) { return std::acos(1.0 / value); }
    static double arccsc(double value) { return std::asin(1.0 / value); }

    // Hyperbolic functions
    static double sinh(double angle) { return std::sinh(angle); }
    static double cosh(double angle) { return std::cosh(angle); }
    static double tanh(double angle) { return std::tanh(angle); }
    static double coth(double angle) { return 1.0 / std::tanh(angle); }
    static double sech(double angle) { return 1.0 / std::cosh(angle); }
    static double csch(double angle) { return 1.0 / std::sinh(angle); }

    // Hyperbolic inverse functions
    static double arcsinh(double value) { return std::asinh(value); }
    static double arccosh(double value) { return std::acosh(value); }
    static double arctanh(double value) { return std::atanh(value); }
    static double arccoth(double value) { return std::atanh(1.0 / value); }
    static double arcsech(double value) { return std::acosh(1.0 / value); }
    static double arccsch(double value) { return std::asinh(1.0 / value); }

    // Degree to radian conversion
    static double toRadians(double degrees) { return degrees * DEG_TO_RAD; }

    // Radian to degree conversion
    static double toDegrees(double radians) { return radians * RAD_TO_DEG; }

    // Calculate the hypotenuse of a right triangle given its legs
    static double hypotenuse(double a, double b) { return std::hypot(a, b); }

    // Calculate the length of an arc given the radius and central angle (in radians)
    static double arcLength(double radius, double angle) { return radius * angle; }

    // Calculate the area of a sector of a circle given the radius and central angle (in radians)
    static double sectorArea(double radius, double angle) { return 0.5 * radius * radius * angle; }

    // Calculate the distance between two points in 2D Cartesian coordinates
    static double distance(double x1, double y1, double x2, double y2) {
        double dx = x2 - x1;
        double dy = y2 - y1;
        return std::hypot(dx, dy);
    }
};

class Advanced_Trigonometry {
public:
    // Constants
    static constexpr double PI = 3.14159265358979323846;
    static constexpr double DEG_TO_RAD = PI / 180.0;
    static constexpr double RAD_TO_DEG = 180.0 / PI;
    static constexpr double EPSILON = std::numeric_limits<double>::epsilon(); // Machine epsilon

    // Taylor series expansions
    // Calculate sine using Taylor series expansion
    static double sinTaylor(double angle) {
        double result = 0.0;
        double term = angle;
        for (int n = 1; std::abs(term) >= EPSILON; ++n) {
            result += term;
            term *= -angle * angle / ((2 * n) * (2 * n + 1));
        }
        return result;
    }

    // Calculate cosine using Taylor series expansion
    static double cosTaylor(double angle) {
        double result = 0.0;
        double term = 1.0;
        for (int n = 0; std::abs(term) >= EPSILON; ++n) {
            result += term;
            term *= -angle * angle / ((2 * n + 2) * (2 * n + 1));
        }
        return result;
    }

    // Calculate tangent using Taylor series expansion (sin/cos)
    static double tanTaylor(double angle) {
        return sinTaylor(angle) / cosTaylor(angle);
    }

    // Calculate cotangent using Taylor series expansion (cos/sin)
    static double cotTaylor(double angle) {
        return cosTaylor(angle) / sinTaylor(angle);
    }

    // Calculate secant using Taylor series expansion (1/cos)
    static double secTaylor(double angle) {
        return 1.0 / cosTaylor(angle);
    }

    // Calculate cosecant using Taylor series expansion (1/sin)
    static double cscTaylor(double angle) {
        return 1.0 / sinTaylor(angle);
    }

    // Calculate the factorial of a non-negative integer
    static double factorial(int n) {
        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }

    // Calculate the permutation of n objects taken r at a time
    static double permutation(int n, int r) {
        return factorial(n) / factorial(n - r);
    }

    // Calculate the combination of n objects taken r at a time
    static double combination(int n, int r) {
        return permutation(n, r) / factorial(r);
    }
};