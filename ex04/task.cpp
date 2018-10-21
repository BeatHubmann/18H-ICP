#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include <tuple>

void PopulateLattice(int lattice[], const int size, const double treshold, const int seed)
{
    std::mt19937 mt_rand(seed);
    std::uniform_real_distribution<double> uniform_zero_one(0.0, 1.0);

    for (auto ij= 0; ij < size*size; ij++)
    {
        if (uniform_zero_one(mt_rand) < treshold) lattice[ij]= 1; 
    }
}

int FindOriginalCluster(const int m_k_cluster_mass[], const int candidate_k)
{
    if (m_k_cluster_mass[candidate_k] > -1) // To intercept endless loop when falsely called with an original cluster
        return candidate_k;
    const int parent_k= -m_k_cluster_mass[candidate_k];
    if (m_k_cluster_mass[parent_k] > 0)
        return parent_k;
    else
        return FindOriginalCluster(m_k_cluster_mass, parent_k);
}

void CheckSite(int lattice[], int m_k_cluster_mass[], int& k_cluster_label,
               const int size, const long int site_index)
{
    const int top_neighbor= (site_index / size > 0 ? lattice[site_index - size] : 0); // not in top row
    const int left_neighbor= (site_index % size > 0 ? lattice[site_index - 1] : 0); // not at left border
    if (top_neighbor == 0 && left_neighbor == 0) // Start a new cluster
    {
        k_cluster_label++;
        lattice[site_index]= k_cluster_label;
        m_k_cluster_mass[k_cluster_label]= 1;
    }
    else if (top_neighbor == 0 && left_neighbor > 0) // Tie to cluster on the left
    {
        const int true_k= (m_k_cluster_mass[left_neighbor] < 0 ?
            FindOriginalCluster(m_k_cluster_mass, left_neighbor) : left_neighbor);
        lattice[site_index]= true_k;
        m_k_cluster_mass[true_k]++;
    }
    else if (top_neighbor > 0 && left_neighbor == 0) // Tie to cluster on top
    {
        const int true_k= (m_k_cluster_mass[top_neighbor] < 0 ?
            FindOriginalCluster(m_k_cluster_mass, top_neighbor) : top_neighbor);
        lattice[site_index]= true_k;
        m_k_cluster_mass[true_k]++;
    }
    else if (top_neighbor > 0 && left_neighbor > 0 && top_neighbor != left_neighbor) // Choose cluster
    {
        const int adopting_cluster= std::min(top_neighbor, left_neighbor);
        const int adopted_cluster= std::max(top_neighbor, left_neighbor);

        const int true_adopting= (m_k_cluster_mass[adopting_cluster] < 0 ?
            FindOriginalCluster(m_k_cluster_mass, adopting_cluster) : adopting_cluster);
        const int true_adopted= (m_k_cluster_mass[adopted_cluster] < 0 ?
            FindOriginalCluster(m_k_cluster_mass, adopted_cluster) : adopted_cluster);
        lattice[site_index]= true_adopting;
        if (true_adopting != true_adopted)
        {
            m_k_cluster_mass[true_adopting] += (m_k_cluster_mass[true_adopted] + 1);
            m_k_cluster_mass[true_adopted]= -true_adopting;
        }
        else if (true_adopting == true_adopted)
            m_k_cluster_mass[true_adopting]++;
    }
    else if (top_neighbor > 0 && left_neighbor > 0 && top_neighbor == left_neighbor) // Tie to top&left cluster
    {
        const int true_k= (m_k_cluster_mass[top_neighbor] < 0 ?
            FindOriginalCluster(m_k_cluster_mass, top_neighbor) : top_neighbor);
        lattice[site_index]= true_k;
        m_k_cluster_mass[true_k]++;
    }
    // To print lattice for quick debugging 
    // for (int i=0; i < size; i++)
    //     {
    //         for (int j=0; j < size; j++)
    //             std::cout << lattice[i*size + j]<< "/" <<m_k_cluster_mass[lattice[i*size + j]] << "\t";
    //         std::cout << "\n";
    //     }
    // std::cout << "\n**********************************************\n";
}

void TagLargestCluster(int lattice[], const int m_k_cluster_mass[], const int size, const int largest_cluster_id)
{
    for (auto ij= 0; ij < size*size; ij++)
    {
        if (lattice[ij] > 0)
        {
            if (lattice[ij] == largest_cluster_id || FindOriginalCluster(m_k_cluster_mass, lattice[ij]) == largest_cluster_id)
                lattice[ij]= 1;
            else
                lattice[ij]= 0;
        }
    }

}

// Returns center index of occupied cells in box if any, else -1 as well
// as total mass of box
const std::tuple<int, int> ScanBox(const int lattice[], const int size, const int center_index, const int side_length)
{
        if (center_index < 0) // Intercept being called with nonvalid center
            return std::make_tuple(-1, 0);
        
        const int search_range{side_length / 2};
        const int center_row{center_index / size};
        const int center_col{center_index % size};
        const int start_row{std::max(0, center_row - search_range)};
        const int start_col{std::max(0, center_col - search_range)};
        const int end_row{std::min(size, start_row + side_length)};
        const int end_col{std::min(size, start_col + side_length)};

        std::vector<int> occupied{-1}; // To store indices of occupied sites with -1 default

        for (auto i= start_row; i < end_row; i++)
            for (auto j= start_col; j < end_col; j++)
                if (lattice[i * size + j] > 0) // current position
                    occupied.push_back(i * size + j);   

        return std::make_tuple(occupied.at(occupied.size() / 2), occupied.size() -1); // box mass is occupied.size() minus -1 default
}


int FindOccupiedNearCenter(const int lattice[], const int size, int search_size)
{
    const int center{(size / 2) * size + size / 2};
    if (lattice[center] == 1)
        return center;

    int center_index{-1}; //  Setup return variable
    while (center_index < 0 && search_size < size)
    {
        center_index = std::get<0>(ScanBox(lattice, size, center, search_size));
        search_size++;
    }

    return center_index;
}

const std::tuple<int, int> SandboxAlgorithm(const int lattice[], const int size, const int center, const int range)
{
    return std::make_tuple(range, std::get<1>(ScanBox(lattice, size, center, range)));
}


const std::tuple<int, int> BoxCountAlgorithm(const int lattice[], const int size, const int box_size)
{
    std::vector<int> box_centers;
    const int start{box_size / 2};
    for (auto i= start; i < size; i+= box_size)
        for (auto j= start; j < size; j+= box_size)
            box_centers.push_back(i * size + j);
    
    int num_non_empty_boxes{0};
    for (auto& box_center : box_centers)
        if (std::get<1>(ScanBox(lattice, size, box_center, box_size)) > 0)
            num_non_empty_boxes++;

    return std::make_tuple(box_size, num_non_empty_boxes);
}

std::vector<std::tuple<int, int>> RunExperiment1(const int size, const double treshold, const int seed)
{
    int* lattice= (int*)std::calloc(size*size, sizeof(int));
    PopulateLattice(lattice, size, treshold, seed);

    int k_cluster_label= 2;
    int* m_k_cluster_mass= (int*)std::calloc(size*size, sizeof(int));

    for (auto ij= 0; ij < size*size; ij++)
        if (lattice[ij] == 1)
            CheckSite(lattice, m_k_cluster_mass, k_cluster_label, size, ij);        

    int largest_cluster_id{0};
    for (auto k= 2; k < size*size; k++) // Look for largest cluster
        if (m_k_cluster_mass[k] > m_k_cluster_mass[largest_cluster_id])
            largest_cluster_id= k;

    TagLargestCluster(lattice, m_k_cluster_mass, size, largest_cluster_id);

    const int search_size{std::max(size / 500, 2)}; // Heuristic estimate
    const int start_box_range{3}; // Given by task
    const int center= FindOccupiedNearCenter(lattice, size, search_size);

    std::vector<std::tuple<int, int>> mass_in_range;
    for (auto range= 3; range <= size; range++)
        mass_in_range.push_back(SandboxAlgorithm(lattice, size, center, range));

    std::free(m_k_cluster_mass);
    std::free(lattice);

    return mass_in_range;
}

std::vector<std::tuple<int, int>> RunExperiment2(const int size, const double treshold, const int seed)
{
    int* lattice= (int*)std::calloc(size*size, sizeof(int));
    PopulateLattice(lattice, size, treshold, seed);

    int k_cluster_label= 2;
    int* m_k_cluster_mass= (int*)std::calloc(size*size, sizeof(int));

    for (auto ij= 0; ij < size*size; ij++)
        if (lattice[ij] == 1)
            CheckSite(lattice, m_k_cluster_mass, k_cluster_label, size, ij);        

    int largest_cluster_id{0};
    for (auto k= 2; k < size*size; k++) // Look for largest cluster
        if (m_k_cluster_mass[k] > m_k_cluster_mass[largest_cluster_id])
            largest_cluster_id= k;

    TagLargestCluster(lattice, m_k_cluster_mass, size, largest_cluster_id);

    std::vector<std::tuple<int, int>> non_empty_boxes;
    for (auto box_size= size; box_size > 0; box_size /= 2)
        non_empty_boxes.push_back(BoxCountAlgorithm(lattice, size, box_size));

    std::free(m_k_cluster_mass);
    std::free(lattice);

    return non_empty_boxes;
}

int main(int argc, char* argv[])
{
    if (argc < 3 || argc > 3 || argv[1] == "-h") // check command line arguments and give some help
     {  
        std::cerr << "Usage: " << argv[0]
                  << " int_lattice_side_length int_initial_seed \n" << std::endl;
        return 1;
     }
    const int size= atoi(argv[1]);
    const int initial_seed= atoi(argv[2]);

    std::vector<double> tresholds{0.58, 0.592746, 0.61}; // Init vector of probs and make sure p_c is included

    for (auto& treshold : tresholds) // Run experiments on all chosen occupation probs
    {
        std::cout << "Experiment_1:Sandbox_Algorithm" << std::endl;

        std::vector<std::tuple<int, int>> results= RunExperiment1(size, treshold, initial_seed);

        std::cout << "p=" << treshold << std::endl;
       
        for (auto& result : results)
            std::cout << std::get<0>(result) << "\t" << std::get<1>(result) << std::endl;
        std::cout << "\n\n";
    }

    for (auto& treshold : tresholds)
    {
        std::cout << "Experiment_2:Box_Counting_Algorithm" << std::endl;

        std::vector<std::tuple<int, int>> results= RunExperiment2(size, treshold, initial_seed);

        std::cout << "p=" << treshold << std::endl;
       
        for (auto& result : results)
            std::cout << std::get<0>(result) << "\t" << std::get<1>(result) << std::endl;
        std::cout << "\n\n";
    }

    return 0;
}
