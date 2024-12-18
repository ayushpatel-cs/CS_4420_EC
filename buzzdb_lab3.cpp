#include <iostream>
#include <algorithm>
#include <fstream>
#include <chrono>
#include <map>
#include <vector>
#include <list>
#include <unordered_map>
#include <string>
#include <memory>
#include <sstream>
#include <limits>
#include <thread>
#include <queue>
#include <optional>
#include <random>
#include <mutex>
#include <shared_mutex>
#include <cassert>
#include <cstring> 
#include <exception>
#include <atomic>
#include <set>

#define UNUSED(p)  ((void)(p))

#define ASSERT_WITH_MESSAGE(condition, message) \
    do { \
        if (!(condition)) { \
            std::cerr << "Assertion \033[1;31mFAILED\033[0m: " << message << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            std::abort(); \
        } \
    } while(0)

enum FieldType { INT, FLOAT, STRING };

std::string applyRLE(const std::string& input) {
    std::stringstream ss;
    for (size_t i = 0; i < input.length(); i++) {
        size_t count = 1;
        while (i + 1 < input.length() && input[i] == input[i + 1]) {
            i++;
            count++;
        }
        ss << input[i] << count;
    }
    return ss.str();
}

std::string reverseRLE(const std::string& input) {
    std::stringstream ss;
    for (size_t i = 0; i < input.length(); i++) {
        char ch = input[i];
        i++;
        std::string count_str;
        while (i < input.length() && std::isdigit(input[i])) {
            count_str += input[i];
            i++;
        }
        i--;
        int count = std::stoi(count_str);
        ss << std::string(count, ch);
    }
    return ss.str();
}

// Define a basic Field variant class that can hold different types
class Field {
public:
    FieldType type;
    std::unique_ptr<char[]> data;
    size_t data_length;
    bool is_compressed = false; 


public:


    Field(int i) : type(INT) { 
        data_length = sizeof(int);
        data = std::make_unique<char[]>(data_length);
        std::memcpy(data.get(), &i, data_length);
    }

    Field(float f) : type(FLOAT) { 
        data_length = sizeof(float);
        data = std::make_unique<char[]>(data_length);
        std::memcpy(data.get(), &f, data_length);
    }

    Field(const std::string& s) : type(STRING) {
        data_length = s.size() + 1;  
        data = std::make_unique<char[]>(data_length);
        std::memcpy(data.get(), s.c_str(), data_length);
    }

    Field& operator=(const Field& other) {
        if (&other == this) {
            return *this;
        }
        type = other.type;
        data_length = other.data_length;
        std::memcpy(data.get(), other.data.get(), data_length);
        return *this;
    }

    Field(Field&& other){
        type = other.type;
        data_length = other.data_length;
        std::memcpy(data.get(), other.data.get(), data_length);
    }

    FieldType getType() const { return type; }
    int asInt() const { 
        return *reinterpret_cast<int*>(data.get());
    }
    float asFloat() const { 
        return *reinterpret_cast<float*>(data.get());
    }
    std::string asString() const { 
        return std::string(data.get());
    }

    std::string serialize() {
        std::stringstream buffer;
        buffer << type << ' ' << data_length << ' ' << is_compressed << ' ';
        if (type == STRING) {
            buffer << data.get() << ' ';
        } else if (type == INT) {
            buffer << *reinterpret_cast<int*>(data.get()) << ' ';
        } else if (type == FLOAT) {
            buffer << *reinterpret_cast<float*>(data.get()) << ' ';
        }
        return buffer.str();
    }

    void serialize(std::ofstream& out) {
        std::string serializedData = this->serialize();
        out << serializedData;
    }

    static std::unique_ptr<Field> deserialize(std::istream& in) {
        int type; in >> type;
        size_t length; in >> length;
        bool is_compressed; in >> is_compressed;

        std::unique_ptr<Field> field;
        if (type == STRING) {
            std::string val; in >> val;
            field = std::make_unique<Field>(val);
        } else if (type == INT) {
            int val; in >> val;
            field = std::make_unique<Field>(val);
        } else if (type == FLOAT) {
            float val; in >> val;
            field = std::make_unique<Field>(val);
        }
        
        if (field) {
            field->is_compressed = is_compressed;
        }
        return field;
    }

    void print() const{
        switch(getType()){
            case INT: std::cout << asInt(); break;
            case FLOAT: std::cout << asFloat(); break;
            case STRING: std::cout << asString(); break;
        }
    }

    void compress() {
        if (type == STRING && !is_compressed) {
            std::string original_str = asString();
            std::string compressed_str = applyRLE(original_str);
            if (compressed_str.size() < original_str.size()) {
                data_length = compressed_str.size() + 1;
                data = std::make_unique<char[]>(data_length);
                std::memcpy(data.get(), compressed_str.c_str(), data_length);
                is_compressed = true;
            }
        }
    }
    
    void decompress() {

        if (type == STRING && is_compressed) {
            std::string compressed_str(data.get());
            std::string original_str = reverseRLE(compressed_str);
            data_length = original_str.size() + 1;
            data = std::make_unique<char[]>(data_length);
            std::memcpy(data.get(), original_str.c_str(), data_length);
            is_compressed = false;
        }
    }

};

class Tuple {
public:
    std::vector<std::unique_ptr<Field>> fields;

    void addField(std::unique_ptr<Field> field) {
        fields.push_back(std::move(field));
    }

    size_t getSize() const {
        size_t size = 0;
        for (const auto& field : fields) {
            size += field->data_length;
        }
        return size;
    }

    std::string serialize() {
        std::stringstream buffer;
        buffer << fields.size() << ' ';
        for (const auto& field : fields) {
            buffer << field->serialize();
        }
        return buffer.str();
    }

    void serialize(std::ofstream& out) {
        std::string serializedData = this->serialize();
        out << serializedData;
    }

    static std::unique_ptr<Tuple> deserialize(std::istream& in) {
        auto tuple = std::make_unique<Tuple>();
        size_t fieldCount; in >> fieldCount;
        for (size_t i = 0; i < fieldCount; ++i) {
            tuple->addField(Field::deserialize(in));
        }
        return tuple;
    }

    void print() const {
        for (const auto& field : fields) {
            field->print();
            std::cout << " ";
        }
        std::cout << "\n";
    }
};

static constexpr size_t PAGE_SIZE = 4096;  // Fixed page size
static constexpr size_t MAX_SLOTS = 512;   // Fixed number of slots
static constexpr size_t MAX_PAGES= 1000;   // Total Number of pages that can be stored
uint16_t INVALID_VALUE = std::numeric_limits<uint16_t>::max(); // Sentinel value

struct Slot {
    bool empty = true;                 // Is the slot empty?    
    uint16_t offset = INVALID_VALUE;    // Offset of the slot within the page
    uint16_t length = INVALID_VALUE;    // Length of the slot
};

class SlottedPage {
public:
    std::unique_ptr<char[]> page_data = std::make_unique<char[]>(PAGE_SIZE);
    size_t metadata_size = sizeof(Slot) * MAX_SLOTS;

    SlottedPage(){
        Slot* slot_array = reinterpret_cast<Slot*>(page_data.get());
        for (size_t slot_itr = 0; slot_itr < MAX_SLOTS; slot_itr++) {
            slot_array[slot_itr].empty = true;
            slot_array[slot_itr].offset = INVALID_VALUE;
            slot_array[slot_itr].length = INVALID_VALUE;
        }
    }

    bool addTuple(std::unique_ptr<Tuple> tuple) {

        for (auto& field : tuple->fields) {
            field->compress();
        }

        auto serializedTuple = tuple->serialize();
        size_t tuple_size = serializedTuple.size();


        size_t slot_itr = 0;
        Slot* slot_array = reinterpret_cast<Slot*>(page_data.get());        
        for (; slot_itr < MAX_SLOTS; slot_itr++) {
            if (slot_array[slot_itr].empty == true and 
                slot_array[slot_itr].length >= tuple_size) {
                break;
            }
        }
        if (slot_itr == MAX_SLOTS){
            return false;
        }

        slot_array[slot_itr].empty = false;
        size_t offset = INVALID_VALUE;
        if (slot_array[slot_itr].offset == INVALID_VALUE){
            if(slot_itr != 0){
                auto prev_slot_offset = slot_array[slot_itr - 1].offset;
                auto prev_slot_length = slot_array[slot_itr - 1].length;
                offset = prev_slot_offset + prev_slot_length;
            }
            else{
                offset = metadata_size;
            }

            slot_array[slot_itr].offset = offset;
        }
        else{
            offset = slot_array[slot_itr].offset;
        }

        if(offset + tuple_size >= PAGE_SIZE){
            slot_array[slot_itr].empty = true;
            slot_array[slot_itr].offset = INVALID_VALUE;
            return false;
        }

        assert(offset != INVALID_VALUE);
        assert(offset >= metadata_size);
        assert(offset + tuple_size < PAGE_SIZE);

        if (slot_array[slot_itr].length == INVALID_VALUE){
            slot_array[slot_itr].length = tuple_size;
        }

        std::memcpy(page_data.get() + offset, 
                    serializedTuple.c_str(), 
                    tuple_size);

        return true;
    }
    void compress() {
        Slot* slot_array = reinterpret_cast<Slot*>(page_data.get());
        for (size_t slot_itr = 0; slot_itr < MAX_SLOTS; slot_itr++) {
            Slot& slot = slot_array[slot_itr];
            if (!slot.empty) {
                char* tuple_data = page_data.get() + slot.offset;
                std::istringstream iss(std::string(tuple_data, slot.length));
                auto tuple = Tuple::deserialize(iss);

                for (auto& field : tuple->fields) {
                    field->compress();
                }

                std::string serialized_tuple = tuple->serialize();
            }
        }
    }

    void decompress() {
        Slot* slot_array = reinterpret_cast<Slot*>(page_data.get());

        std::unique_ptr<char[]> new_page_data = std::make_unique<char[]>(PAGE_SIZE);
        size_t metadata_end = sizeof(Slot) * MAX_SLOTS;
        size_t current_offset = metadata_end;

        for (size_t slot_itr = 0; slot_itr < MAX_SLOTS; slot_itr++) {
            Slot& slot = slot_array[slot_itr];
            if (!slot.empty) {
                char* tuple_data = page_data.get() + slot.offset;
                std::string serialized_tuple(tuple_data, slot.length);
                std::istringstream iss(serialized_tuple);
                auto tuple = Tuple::deserialize(iss);

                for (auto& field : tuple->fields) {
                    field->decompress();
                }

                std::string decompressed_data = tuple->serialize();
                size_t decompressed_length = decompressed_data.size();

                if (current_offset + decompressed_length > PAGE_SIZE) {
                    std::cerr << "Error: Decompressed data exceeds page size.\n";
                    exit(-1);
                }

                std::memcpy(new_page_data.get() + current_offset, decompressed_data.c_str(), decompressed_length);

                slot.offset = current_offset;
                slot.length = decompressed_length;

                current_offset += decompressed_length;
            }
        }

        std::memcpy(page_data.get(), new_page_data.get(), PAGE_SIZE);
    }


    void deleteTuple(size_t index) {
        Slot* slot_array = reinterpret_cast<Slot*>(page_data.get());
        size_t slot_itr = 0;
        for (; slot_itr < MAX_SLOTS; slot_itr++) {
            if(slot_itr == index and
               slot_array[slot_itr].empty == false){
                slot_array[slot_itr].empty = true;
                break;
               }
        }

        //std::this_thread::sleep_for(std::chrono::milliseconds(100));
    }

    void print() const{
        Slot* slot_array = reinterpret_cast<Slot*>(page_data.get());
        for (size_t slot_itr = 0; slot_itr < MAX_SLOTS; slot_itr++) {
            if (slot_array[slot_itr].empty == false){
                assert(slot_array[slot_itr].offset != INVALID_VALUE);
                const char* tuple_data = page_data.get() + slot_array[slot_itr].offset;
                std::istringstream iss(tuple_data);
                auto loadedTuple = Tuple::deserialize(iss);
                std::cout << "Slot " << slot_itr << " : [";
                std::cout << (uint16_t)(slot_array[slot_itr].offset) << "] :: ";
                loadedTuple->print();
            }
        }
        std::cout << "\n";
    }
};

const std::string database_filename = "buzzdb.dat";
std::optional<uint64_t> groot = std::nullopt;

class StorageManager {
public:    
    std::fstream fileStream;
    size_t num_pages = 0;
    std::mutex io_mutex;

public:
    StorageManager(bool truncate_mode = true){
        auto flags =  truncate_mode ? std::ios::in | std::ios::out | std::ios::trunc 
            : std::ios::in | std::ios::out;
        fileStream.open(database_filename, flags);
        if (!fileStream) {
            // If file does not exist, create it
            fileStream.clear(); // Reset the state
            fileStream.open(database_filename, truncate_mode ? (std::ios::out | std::ios::trunc) : std::ios::out);
        }
        fileStream.close(); 
        fileStream.open(database_filename, std::ios::in | std::ios::out); 

        fileStream.seekg(0, std::ios::end);
        num_pages = fileStream.tellg() / PAGE_SIZE;

        if(num_pages == 0){
            extend();
        }

    }

    ~StorageManager() {
        if (fileStream.is_open()) {
            fileStream.close();
        }
    }

    // Read a page from disk
    std::unique_ptr<SlottedPage> load(uint16_t page_id) {
        fileStream.seekg(page_id * PAGE_SIZE, std::ios::beg);
        auto page = std::make_unique<SlottedPage>();
        // Read the content of the file into the page
        if(fileStream.read(page->page_data.get(), PAGE_SIZE)){
            //std::cout << "Page read successfully from file." << std::endl;
        }
        else{
            std::cerr << "Error: Unable to read data from the file. \n";
            exit(-1);
        }
        return page;
    }

    // Write a page to disk
    void flush(uint16_t page_id, const SlottedPage& page) {
        size_t page_offset = page_id * PAGE_SIZE;        

        // Move the write pointer
        fileStream.seekp(page_offset, std::ios::beg);
        fileStream.write(page.page_data.get(), PAGE_SIZE);        
        fileStream.flush();
    }

    // Extend database file by one page
    void extend() {
        // Create a slotted page
        auto empty_slotted_page = std::make_unique<SlottedPage>();

        // Move the write pointer
        fileStream.seekp(0, std::ios::end);

        // Write the page to the file, extending it
        fileStream.write(empty_slotted_page->page_data.get(), PAGE_SIZE);
        fileStream.flush();

        // Update number of pages
        num_pages += 1;
    }

    void extend(uint64_t till_page_id) {
        std::lock_guard<std::mutex>  io_guard(io_mutex); 
        uint64_t write_size = std::max(static_cast<uint64_t>(0), till_page_id + 1 - num_pages) * PAGE_SIZE;
        if(write_size > 0 ) {
            // std::cout << "Extending database file till page id : "<<till_page_id<<" \n";
            char* buffer = new char[write_size];
            std::memset(buffer, 0, write_size);

            fileStream.seekp(0, std::ios::end);
            fileStream.write(buffer, write_size);
            fileStream.flush();
            
            num_pages = till_page_id+1;
        }
    }

};

using PageID = uint16_t;

class Policy {
public:
    virtual bool touch(PageID page_id) = 0;
    virtual PageID evict() = 0;
    virtual ~Policy() = default;
};

void printList(std::string list_name, const std::list<PageID>& myList) {
        std::cout << list_name << " :: ";
        for (const PageID& value : myList) {
            std::cout << value << ' ';
        }
        std::cout << '\n';
}

class LruPolicy : public Policy {
private:
    // List to keep track of the order of use
    std::list<PageID> lruList;

    // Map to find a page's iterator in the list efficiently
    std::unordered_map<PageID, std::list<PageID>::iterator> map;

    size_t cacheSize;

public:

    LruPolicy(size_t cacheSize) : cacheSize(cacheSize) {}

    bool touch(PageID page_id) override {
        //printList("LRU", lruList);

        bool found = false;
        // If page already in the list, remove it
        if (map.find(page_id) != map.end()) {
            found = true;
            lruList.erase(map[page_id]);
            map.erase(page_id);            
        }

        // If cache is full, evict
        if(lruList.size() == cacheSize){
            evict();
        }

        if(lruList.size() < cacheSize){
            // Add the page to the front of the list
            lruList.emplace_front(page_id);
            map[page_id] = lruList.begin();
        }

        return found;
    }

    PageID evict() override {
        // Evict the least recently used page
        PageID evictedPageId = INVALID_VALUE;
        if(lruList.size() != 0){
            evictedPageId = lruList.back();
            map.erase(evictedPageId);
            lruList.pop_back();
        }
        return evictedPageId;
    }

};

constexpr size_t MAX_PAGES_IN_MEMORY = 10;

class BufferManager {
private:
    using PageMap = std::unordered_map<PageID, SlottedPage>;

    StorageManager storage_manager;
    PageMap pageMap;
    std::unique_ptr<Policy> policy;

public:
    BufferManager(bool storage_manager_truncate_mode = true): 
        storage_manager(storage_manager_truncate_mode),
        policy(std::make_unique<LruPolicy>(MAX_PAGES_IN_MEMORY)) {
            storage_manager.extend(MAX_PAGES);
    }
    
    ~BufferManager() {
        for (auto& pair : pageMap) {
            flushPage(pair.first);
        }
    }

    SlottedPage& fix_page(int page_id) {
        auto it = pageMap.find(page_id);
        if (it != pageMap.end()) {
            policy->touch(page_id);
            return pageMap.find(page_id)->second;
        }

        if (pageMap.size() >= MAX_PAGES_IN_MEMORY) {
            auto evictedPageId = policy->evict();
            if(evictedPageId != INVALID_VALUE){
                // std::cout << "Evicting page " << evictedPageId << "\n";
                storage_manager.flush(evictedPageId, 
                                      pageMap[evictedPageId]);
            }
        }

        auto page = storage_manager.load(page_id);
        policy->touch(page_id);
        // std::cout << "Loading page: " << page_id << "\n";
        pageMap[page_id] = std::move(*page);
        return pageMap[page_id];
    }

    void flushPage(int page_id) {
        auto& page = pageMap[page_id];
        storage_manager.flush(page_id, page);
    }

    void extend(){
        storage_manager.extend();
    }
    
    size_t getNumPages(){
        return storage_manager.num_pages;
    }

};

template<typename KeyT, typename ValueT, typename ComparatorT, size_t PageSize>
class BTree {
    public:
        struct Node {
            /// The level in the tree.
            uint16_t level;

            /// The number of children.
            uint16_t count;

            /// TODO: Add additional members as needed

            // Constructor
            Node(uint16_t level, uint16_t count)
                : level(level), count(count) {}

            /// Is the node a leaf node?
            bool is_leaf() const { return level == 0; }
        };

        struct InnerNode: public Node {
            /// The capacity of a node.
            /// TODO think about the capacity that the nodes have.
            static constexpr uint32_t kCapacity = 82;

            /// The keys.
            KeyT keys[kCapacity - 1];

            /// The children.
            uint64_t children[kCapacity];

            /// Constructor.
            InnerNode() : Node(0, 0) {}

            /// Get the index of the first key that is not less than than a provided key.
            /// @param[in] key          The key that should be searched.
            std::pair<uint32_t, bool> lower_bound(const KeyT &key) {
            // TODO: remove the below lines of code 
                // and add your implementation here
                int32_t left = 0;
                int32_t right = this->count - 2; 
                ComparatorT cmp;

                while (left <= right) {
                    uint32_t mid = left + (right - left) / 2;
                    if (cmp(keys[mid], key)) {
                        left = mid + 1;
                    } else if (cmp(key, keys[mid])) {
                        right = mid - 1;
                    } else {
                        return { mid, true };
                    }
                }
                return { left, false };
            }

            uint32_t find_child_index(const KeyT& key) const {
                int32_t left = 0;
                int32_t right = this->count - 2; 
                ComparatorT cmp;

                while (left <= right) {
                    uint32_t mid = left + (right - left) / 2;
                    if (cmp(key, keys[mid])) {
                        right = mid - 1;
                    } else {
                        left = mid + 1;
                    }
                }
                return left;
            }
            /// Insert a key.
            /// @param[in] key          The separator that should be inserted.
            /// @param[in] split_page   The id of the split page that should be inserted.
            void insert(const KeyT &key, uint64_t split_page) {
            // TODO: remove the below lines of code 
                // and add your implementation here
                assert(this->count < kCapacity); 

                auto [index, found] = lower_bound(key);

                for (int32_t i = this->count - 1; i >= static_cast<int32_t>(index); --i) {
                    keys[i + 1] = keys[i];
                }

                for (int32_t i = this->count; i >= static_cast<int32_t>(index + 1); --i) {
                    children[i + 1] = children[i];
                }

                keys[index] = key;
                index += 1;
                children[index] = split_page;

                this->count = this->count + 1;
            
            }

            /// Split the inner node.
            /// @param[in] inner_node       The inner node being split.
            /// @return                 The separator key.
            KeyT split(InnerNode* inner_node) {
                // TODO: remove the below lines of code 
                // and add your implementation here
                uint32_t mid = (this->count + 1) / 2;
                KeyT separator = keys[mid];

                uint32_t right_key_count = this->count - mid - 1;
                for (uint32_t i = 0; i < right_key_count; ++i) {
                    inner_node->keys[i] = keys[mid + 1 + i];
                }

                uint32_t right_child_count = right_key_count + 1;
                for (uint32_t i = 0; i < right_child_count; ++i) {
                    inner_node->children[i] = children[mid + 1 + i];
                }

                inner_node->count = right_child_count;
                inner_node->level = this->level;

                this->count = mid;
                this->count = this->count + 1;

                return separator;
            }

        };

        struct LeafNode: public Node {
            /// The capacity of a node.
            /// TODO think about the capacity that the nodes have.
            static constexpr uint32_t kCapacity =82;

            /// The keys.
            KeyT keys[kCapacity + 1];

            /// The values.
            ValueT values[kCapacity + 1];

            /// Constructor.
            LeafNode() : Node(0, 0) {}

            std::pair<uint32_t, bool> lower_bound(const KeyT& key) const {
                int32_t left = 0;
                int32_t right = this->count - 1;
                ComparatorT cmp;

                while (left <= right) {
                    uint32_t mid = left + (right - left) / 2;
                    if (cmp(keys[mid], key)) {
                        left = mid + 1;
                    } else if (cmp(key, keys[mid])) {
                        right = mid - 1;
                    } else {
                        return { mid, true };
                    }
                }
                return { left, false };
            }

            /// Insert a key.
            /// @param[in] key          The key that should be inserted.
            /// @param[in] value        The value that should be inserted.
            void insert(const KeyT &key, const ValueT &value) {
                // HINT
                auto [index, found] = lower_bound(key);
                assert(this->count < kCapacity + 1); 
                if (found) {
                    values[index] = value;
                } else {
                    for (int32_t i = this->count - 1; i >= static_cast<int32_t>(index); --i) {
                assert(static_cast<uint32_t>(i) + 1 < kCapacity + 1);
                        keys[i + 1] = keys[i];
                        values[i + 1] = values[i];
                    }
                    keys[index] = key;
                    values[index] = value;
                    this->count = this->count + 1;
                }
                assert(this->count <= kCapacity + 1);
            }

            /// Erase a key.
            void erase(const KeyT &key) {
                auto [index, found] = lower_bound(key);
                if (found) {
                    for (uint32_t i = index; i < (uint64_t) (this->count - 1); ++i) {
                        keys[i] = keys[i + 1];
                        values[i] = values[i + 1];
                    }
                    this->count = this->count -1;
                }
            }

            /// Split the leaf node.
            /// @param[in] leaf_node       The leaf node being split
            /// @return                 The separator key.
            KeyT split(LeafNode* leaf_node) {
                
                uint32_t mid = (this->count + 1) / 2;
                uint32_t num_keys_to_move = this->count - mid;

                assert(mid <= kCapacity);
                assert(num_keys_to_move <= kCapacity);

                for (uint32_t i = 0; i < num_keys_to_move; ++i) {
                    leaf_node->keys[i] = keys[mid + i];
                    leaf_node->values[i] = values[mid + i];
                }

                leaf_node->count = num_keys_to_move;
                leaf_node->level = this->level;

                this->count = mid;


                return leaf_node->keys[0];
            }

        };

        /// The root.
        std::optional<uint64_t> root;

        /// The buffer manager
        BufferManager& buffer_manager;

        /// Next page id.
        /// You don't need to worry about about the page allocation.
        /// (Neither fragmentation, nor persisting free-space bitmaps)
        /// Just increment the next_page_id whenever you need a new page.
        uint64_t next_page_id;

        /// Constructor.
        BTree(BufferManager &buffer_manager): buffer_manager(buffer_manager) {
            // TODO
            if (buffer_manager.getNumPages() == 0) {
            buffer_manager.extend();
            }
            root = std::nullopt;
            next_page_id = 1;
            if (groot) {
                next_page_id = static_cast<uint64_t>(*groot) + 1;
            }
            auto& meta_page = buffer_manager.fix_page(0);
            Slot* slot_array = reinterpret_cast<Slot*>(meta_page.page_data.get());

            for (size_t slot_itr = 0; slot_itr < MAX_SLOTS; slot_itr++) {

                if (!slot_array[slot_itr].empty) {
                    char* tuple_data = meta_page.page_data.get() + slot_array[slot_itr].offset;
                    std::istringstream iss(std::string(tuple_data, slot_array[slot_itr].length));
                    auto tuple = Tuple::deserialize(iss);
                    break;
                }
            }
        }

        // ~BTree() {
        //     if (root.has_value()) {
        //         auto& meta_page = buffer_manager.fix_page(0);
        //         std::memset(meta_page.page_data.get(), 0, PAGE_SIZE);
        //         auto root_tuple = std::make_unique<Tuple>();
        //         root_tuple->addField(std::make_unique<Field>(static_cast<int>(*root)));
        //         SlottedPage slotted_page;
        //         buffer_manager.flushPage(0);
        //     }
        // }

        /// Lookup an entry in the tree.
        /// @param[in] key      The key that should be searched.
        std::optional<ValueT> lookup(const KeyT &key) {
            root = groot;
            if (!root.has_value()) return std::nullopt;
            return lookup_rec(*root, key);
        }

        std::optional<ValueT> lookup_rec(uint64_t node_page_id, const KeyT& key) {

            // std::cout << "Test1" << std::endl;
            // std::cout << node_page_id << " " << key << std::endl;


            auto& node_page = buffer_manager.fix_page(node_page_id);
            auto node = reinterpret_cast<Node*>(node_page.page_data.get());
            // std::cout << "Test2" << std::endl;

            if (node->is_leaf()) {
                auto leaf_node = reinterpret_cast<LeafNode*>(node);
                auto [index, found] = leaf_node->lower_bound(key);
                // std::cout << "Test3A" << std::endl;
                if (found) {
                    return leaf_node->values[index];
                } else {
                    return std::nullopt;
                }
            } else {
                auto inner_node = reinterpret_cast<InnerNode*>(node);
                uint32_t index = inner_node->find_child_index(key);
                uint64_t child_page_id = inner_node->children[index];

                return lookup_rec(child_page_id, key);
            }

        }

        /// Erase an entry in the tree.
        /// @param[in] key      The key that should be searched.
        void erase(const KeyT &key) {
                // TODO
            if (!root.has_value()) {
                return;
            }
            erase_rec(*root, key);
            auto& root_page = buffer_manager.fix_page(*root);
            auto root_node = reinterpret_cast<Node*>(root_page.page_data.get());
            if (!root_node->is_leaf()) {
                if (root_node->count ==0){
                    root = reinterpret_cast<InnerNode*>(root_node)->children[0];
                    groot = root;
                }
            }
        }

        void erase_rec(uint64_t node_page_id, const KeyT& key) {
            auto& node_page = buffer_manager.fix_page(node_page_id);
            auto node = reinterpret_cast<Node*>(node_page.page_data.get());

            if (node->is_leaf()) {
                auto leaf_node = reinterpret_cast<LeafNode*>(node);
                leaf_node->erase(key);
                buffer_manager.flushPage(node_page_id);
            } else {
                auto inner_node = reinterpret_cast<InnerNode*>(node);
                uint32_t index = inner_node->find_child_index(key);
                uint64_t child_page_id = inner_node->children[index];
                erase_rec(child_page_id, key);
                buffer_manager.flushPage(node_page_id);
            }
        }

        /// Inserts a new entry into the tree.
        /// @param[in] key      The key that should be inserted.
        /// @param[in] value    The value that should be inserted.
        void insert(const KeyT &key, const ValueT &value) {
             if (!root.has_value()) {
                root = next_page_id++;
                
                auto& root_page = buffer_manager.fix_page(*root);
                auto rn = reinterpret_cast<LeafNode*>(root_page.page_data.get());
                new (rn) LeafNode();
                rn->insert(key, value);
                buffer_manager.flushPage(*root);
            } else {
                auto nc = insert_recursive(*root, key, value);
                if (nc.has_value()) {
                    auto old_root_page_id = *root;
                    root = next_page_id++;
                    auto& root_page = buffer_manager.fix_page(*root);
                    auto rn = reinterpret_cast<InnerNode*>(root_page.page_data.get());
                    new (rn) InnerNode();
                    rn->level = reinterpret_cast<Node*>(buffer_manager.fix_page(old_root_page_id).page_data.get())->level + 1;
                    rn->keys[0] = nc->first;
                    rn->children[0] = old_root_page_id;
                    rn->children[1] = nc->second;
                    rn->count = 2;
                    buffer_manager.flushPage(*root);
                }
            }
            groot = root;
        }

        std::optional<std::pair<KeyT, uint64_t>> insert_recursive(uint64_t node_page_id, const KeyT& key, const ValueT& value) {
            auto& node_page = buffer_manager.fix_page(node_page_id);
            auto node = reinterpret_cast<Node*>(node_page.page_data.get());

            if (node->is_leaf()) {
                auto leaf_node = reinterpret_cast<LeafNode*>(node);

                if (leaf_node->count >= LeafNode::kCapacity + 1) {
                    std::abort();
                }
                if (leaf_node->count >= LeafNode::kCapacity) {

                    auto new_child_page_id = next_page_id++;
                    auto& new_leaf_page = buffer_manager.fix_page(new_child_page_id);
                    auto new_leaf_node = reinterpret_cast<LeafNode*>(new_leaf_page.page_data.get());
                    new (new_leaf_node) LeafNode();
                    KeyT separator = leaf_node->split(new_leaf_node);


                    if (!ComparatorT()(key, separator)) {
                        new_leaf_node->insert(key, value);
                    } else {
                        leaf_node->insert(key, value);
                    }

                    buffer_manager.flushPage(node_page_id);
                    buffer_manager.flushPage(new_child_page_id);

                    return std::make_pair(separator, new_child_page_id);
                } else {
                    leaf_node->insert(key, value);
                    buffer_manager.flushPage(node_page_id);
                    return std::nullopt;
                }
            } else {
                auto inner_node = reinterpret_cast<InnerNode*>(node);
                uint32_t index = inner_node->find_child_index(key);
                uint64_t child_page_id = inner_node->children[index];

                auto new_child = insert_recursive(child_page_id, key, value);
                if (new_child.has_value()) {
                    if (inner_node->count >= InnerNode::kCapacity) {

                        auto new_inner_page_id = next_page_id++;
                        auto& new_inner_page = buffer_manager.fix_page(new_inner_page_id);
                        auto new_inner_node = reinterpret_cast<InnerNode*>(new_inner_page.page_data.get());
                        new (new_inner_node) InnerNode();
                        KeyT separator = inner_node->split(new_inner_node);
                        if (ComparatorT()(new_child->first, separator)) {
                            inner_node->insert(new_child->first, new_child->second);
                        } else {
                            new_inner_node->insert(new_child->first, new_child->second);
                        }

                        buffer_manager.flushPage(node_page_id);
                        buffer_manager.flushPage(new_inner_page_id);

                        return std::make_pair(separator, new_inner_page_id);
                    } else {
                        inner_node->insert(new_child->first, new_child->second);
                        buffer_manager.flushPage(node_page_id);
                        return std::nullopt;
                    }
                } else {
                    buffer_manager.flushPage(node_page_id);
                    return std::nullopt;
                }
            }
        }
};

int main(int argc, char* argv[]) {
    auto testApplyRLE = []() {
        std::string input = "aaabbbcccaaa";
        std::string expected = "a3b3c3a3";
        std::string output = applyRLE(input);
        assert(output == expected);
        std::cout << "testApplyRLE passed.\n";
    };

    auto testReverseRLE = []() {
        std::string input = "a3b3c3a3";
        std::string expected = "aaabbbcccaaa";
        std::string output = reverseRLE(input);
        assert(output == expected);
        std::cout << "testReverseRLE passed.\n";
    };

    auto testRLEEdgeCases = []() {
        std::string input = "";
        std::string output = applyRLE(input);
        assert(output == "");
        output = reverseRLE(output);
        assert(output == "");
        std::cout << "testRLEEdgeCases (empty string) passed.\n";

        // Test string with no repeats
        input = "abcdef";
        output = applyRLE(input);
        // Each character should have count 1
        assert(output == "a1b1c1d1e1f1");
        output = reverseRLE(output);
        assert(output == "abcdef");
        std::cout << "testRLEEdgeCases (no repeats) passed.\n";
    };

    auto testFieldCompression = []() {
    std::string original_str = "aaabbbcccaaa";
    Field string_field(original_str);

    string_field.compress();
    assert(string_field.is_compressed == true);

    std::string compressed_str = string_field.asString();
    assert(compressed_str == "a3b3c3a3");

    string_field.decompress();
    assert(string_field.is_compressed == false);
    std::string decompressed_str = string_field.asString();
    assert(decompressed_str == original_str);

    std::cout << "STRING compression test passed.\n";

    int original_int = 42;
    Field int_field(original_int);

    int_field.compress(); 
    assert(int_field.asInt() == original_int); 
    assert(int_field.is_compressed == false); 

    int_field.decompress(); 
    assert(int_field.asInt() == original_int);

    std::cout << "INT compression test passed.\n";

    float original_float = 3.14159f;
    Field float_field(original_float);

    float_field.compress(); 
    assert(float_field.asFloat() == original_float); 
    assert(float_field.is_compressed == false); 

    float_field.decompress(); 
    assert(float_field.asFloat() == original_float);

    std::cout << "FLOAT compression test passed.\n";
};

    auto testTupleInsertionAndRetrieval = []() {
        SlottedPage page;

        auto tuple = std::make_unique<Tuple>();
        tuple->addField(std::make_unique<Field>(42)); 
        tuple->addField(std::make_unique<Field>("aaabbbcccaaa")); 

        bool success = page.addTuple(std::move(tuple));
        assert(success == true);

        Slot* slot_array = reinterpret_cast<Slot*>(page.page_data.get());
        for (size_t slot_itr = 0; slot_itr < MAX_SLOTS; slot_itr++) {
            if (!slot_array[slot_itr].empty) {
                char* tuple_data = page.page_data.get() + slot_array[slot_itr].offset;
                size_t tuple_length = slot_array[slot_itr].length;
                std::istringstream iss(std::string(tuple_data, tuple_length));
                auto loadedTuple = Tuple::deserialize(iss);

                assert(loadedTuple->fields[0]->asInt() == 42);

                loadedTuple->fields[1]->decompress();
                assert(loadedTuple->fields[1]->asString() == "aaabbbcccaaa");
            }
        }

        std::cout << "testTupleInsertionAndRetrieval passed.\n";
    };

    auto testDatabaseOperations = []() {
        BufferManager buffer_manager(true); 
        SlottedPage page;

        for (uint64_t i = 0; i < 10; ++i) {
            std::string str = "aaaabbbbccccddddeeee"; 
            auto tuple = std::make_unique<Tuple>();
            tuple->addField(std::make_unique<Field>(static_cast<int>(i)));
            tuple->addField(std::make_unique<Field>(str)); 

            bool success = page.addTuple(std::move(tuple));
            assert(success == true);
        }

        Slot* slot_array = reinterpret_cast<Slot*>(page.page_data.get());
        for (size_t slot_itr = 0; slot_itr < MAX_SLOTS; slot_itr++) {
            if (!slot_array[slot_itr].empty) {
                char* tuple_data = page.page_data.get() + slot_array[slot_itr].offset;
                std::string serialized_tuple(tuple_data, slot_array[slot_itr].length);
                std::istringstream iss(serialized_tuple);
                auto loadedTuple = Tuple::deserialize(iss);

                int expected_int = static_cast<int>(slot_itr);
                assert(loadedTuple->fields[0]->asInt() == expected_int);

                loadedTuple->fields[1]->decompress();
                std::string expected_str = "aaaabbbbccccddddeeee";
                assert(loadedTuple->fields[1]->asString() == expected_str);
            }
        }

        std::cout << "testDatabaseOperations passed.\n";
    };


    bool execute_all = false;
    std::string selected_test = "-1";

    if(argc < 2) {
        execute_all = true;
    } else {
        selected_test = argv[1];
    }

    using BTree = BTree<uint64_t, uint64_t, std::less<uint64_t>, 1024>;

    // Test 1: InsertEmptyTree
    if(execute_all || selected_test == "1") {
        std::cout<<"...Starting Test 1"<<std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);

        ASSERT_WITH_MESSAGE(tree.root.has_value() == false, 
            "tree.root is not nullptr");

        tree.insert(42, 21);

        ASSERT_WITH_MESSAGE(tree.root.has_value(), 
            "tree.root is still nullptr after insertion");
        
        std::string test = "inserting an element into an empty B-Tree";

        // Fix root page and obtain root node pointer
        SlottedPage* root_page = &buffer_manager.fix_page(*tree.root);
        auto root_node = reinterpret_cast<BTree::Node*>(root_page->page_data.get());

        ASSERT_WITH_MESSAGE(root_node->is_leaf() == true, 
            test + " does not create a leaf node.");
        ASSERT_WITH_MESSAGE(root_node->count == 1,
            test + " does not create a leaf node with count = 1.");

        std::cout << "\033[1m\033[32mPassed: Test 1\033[0m" << std::endl;
    }

    // Test 2: InsertLeafNode
    if(execute_all || selected_test == "2") {
        std::cout<<"...Starting Test 2"<<std::endl;
        BufferManager buffer_manager;
        BTree tree( buffer_manager);
        
        ASSERT_WITH_MESSAGE(tree.root.has_value() == false, 
            "tree.root is not nullptr");

        for (auto i = 0ul; i < BTree::LeafNode::kCapacity; ++i) {
            tree.insert(i, 2 * i);
        }
        ASSERT_WITH_MESSAGE(tree.root.has_value(), 
            "tree.root is still nullptr after insertion");

        std::string test = "inserting BTree::LeafNode::kCapacity elements into an empty B-Tree";

        SlottedPage* root_page = &buffer_manager.fix_page(*tree.root);
        auto root_node = reinterpret_cast<BTree::Node*>(root_page->page_data.get());
        auto root_inner_node = static_cast<BTree::InnerNode*>(root_node);

        ASSERT_WITH_MESSAGE(root_node->is_leaf() == true, 
            test + " creates an inner node as root.");
        ASSERT_WITH_MESSAGE(root_inner_node->count == BTree::LeafNode::kCapacity,
            test + " does not store all elements.");

        std::cout << "\033[1m\033[32mPassed: Test 2\033[0m" << std::endl;
    }

    // Test 3: InsertLeafNodeSplit
    if (execute_all || selected_test == "3") {
        std::cout << "...Starting Test 3" << std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);

        ASSERT_WITH_MESSAGE(tree.root.has_value() == false, 
            "tree.root is not nullptr");

        for (auto i = 0ul; i < BTree::LeafNode::kCapacity; ++i) {
            tree.insert(i, 2 * i);
        }
        ASSERT_WITH_MESSAGE(tree.root.has_value(), 
            "tree.root is still nullptr after insertion");
        
        SlottedPage* root_page = &buffer_manager.fix_page(*tree.root);
        auto root_node = reinterpret_cast<BTree::Node*>(root_page->page_data.get());
        auto root_inner_node = static_cast<BTree::InnerNode*>(root_node);
        
        assert(root_inner_node->is_leaf());
        assert(root_inner_node->count == BTree::LeafNode::kCapacity);

        // Let there be a split...
        tree.insert(424242, 42);

        std::string test =
            "inserting BTree::LeafNode::kCapacity + 1 elements into an empty B-Tree";
        
        ASSERT_WITH_MESSAGE(tree.root.has_value() != false, test + " removes the root :-O");

        SlottedPage* root_page1 = &buffer_manager.fix_page(*tree.root);
        root_node = reinterpret_cast<BTree::Node*>(root_page1->page_data.get());
        root_inner_node = static_cast<BTree::InnerNode*>(root_node);

        ASSERT_WITH_MESSAGE(root_inner_node->is_leaf() == false, 
            test + " does not create a root inner node");
        ASSERT_WITH_MESSAGE(root_inner_node->count == 2, 
            test + " creates a new root with count != 2");

        std::cout << "\033[1m\033[32mPassed: Test 3\033[0m" << std::endl;
    }

    // Test 4: LookupEmptyTree
    if (execute_all || selected_test == "4") {
        std::cout << "...Starting Test 4" << std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);

        std::string test = "searching for a non-existing element in an empty B-Tree";

        ASSERT_WITH_MESSAGE(tree.lookup(42).has_value() == false, 
            test + " seems to return something :-O");

        std::cout << "\033[1m\033[32mPassed: Test 4\033[0m" << std::endl;
    }

    // Test 5: LookupSingleLeaf
    if (execute_all || selected_test == "5") {
        std::cout << "...Starting Test 5" << std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);

        // Fill one page
        for (auto i = 0ul; i < BTree::LeafNode::kCapacity; ++i) {
            tree.insert(i, 2 * i);
            ASSERT_WITH_MESSAGE(tree.lookup(i).has_value(), 
                "searching for the just inserted key k=" + std::to_string(i) + " yields nothing");
        }

        // Lookup all values
        for (auto i = 0ul; i < BTree::LeafNode::kCapacity; ++i) {
            auto v = tree.lookup(i);
            ASSERT_WITH_MESSAGE(v.has_value(), "key=" + std::to_string(i) + " is missing");
            ASSERT_WITH_MESSAGE(*v == 2 * i, "key=" + std::to_string(i) + " should have the value v=" + std::to_string(2 * i));
        }

        std::cout << "\033[1m\033[32mPassed: Test 5\033[0m" << std::endl;
    }

    // Test 6: LookupSingleSplit
    if (execute_all || selected_test == "6") {
        std::cout << "...Starting Test 6" << std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);

        // Insert values
        for (auto i = 0ul; i < BTree::LeafNode::kCapacity; ++i) {
            tree.insert(i, 2 * i);
        }

        tree.insert(BTree::LeafNode::kCapacity, 2 * BTree::LeafNode::kCapacity);
        ASSERT_WITH_MESSAGE(tree.lookup(BTree::LeafNode::kCapacity).has_value(),
            "searching for the just inserted key k=" + std::to_string(BTree::LeafNode::kCapacity + 1) + " yields nothing");

        // Lookup all values
        for (auto i = 0ul; i < BTree::LeafNode::kCapacity + 1; ++i) {
            auto v = tree.lookup(i);
            ASSERT_WITH_MESSAGE(v.has_value(), "key=" + std::to_string(i) + " is missing");
            ASSERT_WITH_MESSAGE(*v == 2 * i, 
                "key=" + std::to_string(i) + " should have the value v=" + std::to_string(2 * i));
        }

        std::cout << "\033[1m\033[32mPassed: Test 6\033[0m" << std::endl;
    }

    // Test 7: LookupMultipleSplitsIncreasing
    if (execute_all || selected_test == "7") {
        std::cout << "...Starting Test 7" << std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);
        auto n = 40 * BTree::LeafNode::kCapacity;

        // Insert values
        for (auto i = 0ul; i < n; ++i) {
            tree.insert(i, 2 * i);
            ASSERT_WITH_MESSAGE(tree.lookup(i).has_value(),
                "searching for the just inserted key k=" + std::to_string(i) + " yields nothing");
        }

        // Lookup all values
        for (auto i = 0ul; i < n; ++i) {
            auto v = tree.lookup(i);
            ASSERT_WITH_MESSAGE(v.has_value(), "key=" + std::to_string(i) + " is missing");
            ASSERT_WITH_MESSAGE(*v == 2 * i,
                "key=" + std::to_string(i) + " should have the value v=" + std::to_string(2 * i));
        }
        std::cout << "\033[1m\033[32mPassed: Test 7\033[0m" << std::endl;
    }

    // Test 8: LookupMultipleSplitsDecreasing
    if (execute_all || selected_test == "8") {
        std::cout << "...Starting Test 8" << std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);
        auto n = 10 * BTree::LeafNode::kCapacity;

        // Insert values
        for (auto i = n; i > 0; --i) {
            tree.insert(i, 2 * i);
            ASSERT_WITH_MESSAGE(tree.lookup(i).has_value(),
                "searching for the just inserted key k=" + std::to_string(i) + " yields nothing");
        }

        // Lookup all values
        for (auto i = n; i > 0; --i) {
            auto v = tree.lookup(i);
            ASSERT_WITH_MESSAGE(v.has_value(), "key=" + std::to_string(i) + " is missing");
            ASSERT_WITH_MESSAGE(*v == 2 * i,
                "key=" + std::to_string(i) + " should have the value v=" + std::to_string(2 * i));
        }

        std::cout << "\033[1m\033[32mPassed: Test 8\033[0m" << std::endl;
    }

    // Test 9: LookupRandomNonRepeating
    if (execute_all || selected_test == "9") {
        std::cout << "...Starting Test 9" << std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);
        auto n = 10 * BTree::LeafNode::kCapacity;

        // Generate random non-repeating key sequence
        std::vector<uint64_t> keys(n);
        std::iota(keys.begin(), keys.end(), n);
        std::mt19937_64 engine(0);
        std::shuffle(keys.begin(), keys.end(), engine);

        // Insert values
        for (auto i = 0ul; i < n; ++i) {
            tree.insert(keys[i], 2 * keys[i]);
            ASSERT_WITH_MESSAGE(tree.lookup(keys[i]).has_value(),
                "searching for the just inserted key k=" + std::to_string(keys[i]) +
                " after i=" + std::to_string(i) + " inserts yields nothing");
        }

        // Lookup all values
        for (auto i = 0ul; i < n; ++i) {
            auto v = tree.lookup(keys[i]);
            ASSERT_WITH_MESSAGE(v.has_value(), "key=" + std::to_string(keys[i]) + " is missing");
            ASSERT_WITH_MESSAGE(*v == 2 * keys[i],
                "key=" + std::to_string(keys[i]) + " should have the value v=" + std::to_string(2 * keys[i]));
        }

        std::cout << "\033[1m\033[32mPassed: Test 9\033[0m" << std::endl;
    }

    // Test 10: LookupRandomRepeating
    if (execute_all || selected_test == "10") {
        std::cout << "...Starting Test 10" << std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);
        auto n = 10 * BTree::LeafNode::kCapacity;

        // Insert & updated 100 keys at random
        std::mt19937_64 engine{0};
        std::uniform_int_distribution<uint64_t> key_distr(0, 99);
        std::vector<uint64_t> values(100);

        for (auto i = 1ul; i < n; ++i) {
            uint64_t rand_key = key_distr(engine);
            values[rand_key] = i;
            tree.insert(rand_key, i);

            auto v = tree.lookup(rand_key);
            ASSERT_WITH_MESSAGE(v.has_value(),
                "searching for the just inserted key k=" + std::to_string(rand_key) +
                " after i=" + std::to_string(i - 1) + " inserts yields nothing");
            ASSERT_WITH_MESSAGE(*v == i,
                "overwriting k=" + std::to_string(rand_key) + " with value v=" + std::to_string(i) +
                " failed");
        }

        // Lookup all values
        for (auto i = 0ul; i < 100; ++i) {
            if (values[i] == 0) {
                continue;
            }
            auto v = tree.lookup(i);
            ASSERT_WITH_MESSAGE(v.has_value(), "key=" + std::to_string(i) + " is missing");
            ASSERT_WITH_MESSAGE(*v == values[i],
                "key=" + std::to_string(i) + " should have the value v=" + std::to_string(values[i]));
        }

        std::cout << "\033[1m\033[32mPassed: Test 10\033[0m" << std::endl;
    }

    // Test 11: Erase
    if (execute_all || selected_test == "11") {
        std::cout << "...Starting Test 11" << std::endl;
        BufferManager buffer_manager;
        BTree tree(buffer_manager);

        // Insert values
        for (auto i = 0ul; i < 2 * BTree::LeafNode::kCapacity; ++i) {
            tree.insert(i, 2 * i);
        }

        // Iteratively erase all values
        for (auto i = 0ul; i < 2 * BTree::LeafNode::kCapacity; ++i) {
            ASSERT_WITH_MESSAGE(tree.lookup(i).has_value(), "k=" + std::to_string(i) + " was not in the tree");
            tree.erase(i);
            ASSERT_WITH_MESSAGE(!tree.lookup(i), "k=" + std::to_string(i) + " was not removed from the tree");
        }
        std::cout << "\033[1m\033[32mPassed: Test 11\033[0m" << std::endl;
    }

    // Test 12: Persistant Btree
    if (execute_all || selected_test == "12") {
        std::cout << "...Starting Test 12" << std::endl;
        unsigned long n =  10 * BTree::LeafNode::kCapacity;

        // Build a tree
        {
            BufferManager buffer_manager;
            BTree tree(buffer_manager);
            
            // Insert values
            for (auto i = 0ul; i < n; ++i) {
                tree.insert(i, 2 * i);
                ASSERT_WITH_MESSAGE(tree.lookup(i).has_value(),
                    "searching for the just inserted key k=" + std::to_string(i) + " yields nothing");
            }

            // Lookup all values
            for (auto i = 0ul; i < n; ++i) {
                auto v = tree.lookup(i);
                ASSERT_WITH_MESSAGE(v.has_value(), "key=" + std::to_string(i) + " is missing");
                ASSERT_WITH_MESSAGE(*v == 2 * i,
                    "key=" + std::to_string(i) + " should have the value v=" + std::to_string(2 * i));
            }
        }

        // recreate the buffer manager and check for existence of the tree
        {
            BufferManager buffer_manager(false);
            BTree tree(buffer_manager);

            // Lookup all values
            for (auto i = 0ul; i < n; ++i) {
                auto v = tree.lookup(i);
                ASSERT_WITH_MESSAGE(v.has_value(), "key=" + std::to_string(i) + " is missing");
                ASSERT_WITH_MESSAGE(*v == 2 * i,
                    "key=" + std::to_string(i) + " should have the value v=" + std::to_string(2 * i));
            }
        }

        std::cout << "\033[1m\033[32mPassed: Test 12\033[0m" << std::endl;



        
    }
    testApplyRLE();
    testReverseRLE();
    testRLEEdgeCases();
    testFieldCompression();
    testTupleInsertionAndRetrieval();
    testDatabaseOperations();

    std::cout << "All tests passed successfully.\n";

    return 0;
}