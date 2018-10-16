#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>
#include <tuple>

void PopulateLattice(int lattice[], const int size, const double treshold, const int seed)
{
    std::mt19937 mt_rand(seed);
    std::uniform_real_distribution<double> uniform_zero_one(0.0, 1.0);

    for (int ij= 0; ij < size*size; ij++)
    {
        if (uniform_zero_one(mt_rand) < treshold) lattice[ij]= 1; 
    }
}

int FindOriginalCluster(int m_k_cluster_mass[], const int candidate_k)
{
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
        else //if (true_adopting == true_adopted)
            m_k_cluster_mass[true_adopting]++;
    }
    else if (top_neighbor > 0 && left_neighbor > 0 && top_neighbor == left_neighbor) // Tie to top&left cluster
    {
        const int true_k= (m_k_cluster_mass[top_neighbor] < 0 ?
            FindOriginalCluster(m_k_cluster_mass, top_neighbor) : top_neighbor);
        lattice[site_index]= true_k;
        m_k_cluster_mass[true_k]++;
    }
    // // To print lattice for quick debugging 
    // for (int i=0; i < size; i++)
    //     {
    //         for (int j=0; j < size; j++)
    //             std::cout << lattice[i*size + j]<< "/" <<m_k_cluster_mass[lattice[i*size + j]] << "\t";
    //         std::cout << "\n";
    //     }
    // std::cout << "\n**********************************************\n";
}

void TagLargestCluster(int lattice[], const int size, const int largest_cluster_id)
{
    for (int ij= 0; ij < size*size; ij++)
    {
        if (lattice[ij] == largest_cluster_id)
            lattice[ij]= 1;
        else
            lattice[ij]= 0;
    }

}

void ScanBox(const int lattice[], const int size, const int center, const int range)
{

    
}


const int FindOccupiedNearCenter(const int lattice[], const int size)
{
    const int starting_index{(size / 2) * size + (size / 2)};
    if (lattice[starting_index] == 1)
        return starting_index;
    int center_index{-1}; // Set to be outside lattice
    int search_range{1};
    while (center_index < 0 && search_range < size / 2)
    {
        const int center_row{center_index % size};
        const int center_col{center_index / size};
        const int start_row{std::max(0, center_row - search_range)};
        const int start_col{std::max(0, center_col - search_range)};
        const int end_row{std::min(size, center_row + search_range + 1)}; // 1 index post the end
        const int end_col{std::min(size, center_col + search_range + 1)}; // 1 index post the end
        int j{start_col};
        while (center_index < 0 && j < e)
        {
            while (center_index < 0 && ij < starting_index + search_range)
            {
                if (lattice[ij] == 1)
                    center_index= ij;
                ij++;
            }
            ij= starting_index - search_range * size - search_range;
        }
         search_range++;
    }
    return center_index;
}

const int SandboxAlgorithm(const int lattice[], const int size, const int center, const int range)
{
    const int starting_row= std::max(0, center % size - size * range/2);
    const int ending_row= std::min(size - 1, center % size + size * range/2);
       for (int i= center - range/2; i <= center)
    return 1;
}

std::vector<std::tuple<int, int>> RunExperiment(const int size, const double treshold, const int seed)
{
    int* lattice= (int*)std::calloc(size*size, sizeof(int));
    PopulateLattice(lattice, size, treshold, seed);

    int k_cluster_label= 2;
    int* m_k_cluster_mass= (int*)std::calloc(size*size, sizeof(int));

    for (long int ij= 0; ij < size*size; ij++)
        if (lattice[ij] == 1)
            CheckSite(lattice, m_k_cluster_mass, k_cluster_label, size, ij);        

    int largest_cluster_id{0};
    int largest_cluster_mass{0};
    for (int k= 2; k < size*size; k++) // Look for largest cluster
        if (m_k_cluster_mass[k] > largest_cluster_mass)
            largest_cluster_id= k;

    TagLargestCluster(lattice, size, largest_cluster_id);

    const int center= FindOccupiedNearCenter(lattice, size);

    std::vector<int> mass_in_radius;

    for (auto radius= 3; radius <= size; radius++)
        mass_in_radius.push_back(SandboxAlgorithm(lattice, size, center, radius));

    std::free(m_k_cluster_mass);
    std::free(lattice);

    return mass_in_radius;
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

    std::vector<double> tresholds{0.57, 0.58, 0.592746, 0.6, 0.61}; // Init vector of probs and make sure p_c is included

    for (auto& treshold : tresholds) // Run experiments on all chosen occupation probs
    {
        std::cout << "#p=" << treshold << std::endl;

        std::vector<int> results= RunExperiment(size, treshold, initial_seed);
       
        for (auto& result : results)
            std::cout << 0 << "\t" << 1 << std::endl;
        std::cout << "\n\n";
    }

    return 0;
}
