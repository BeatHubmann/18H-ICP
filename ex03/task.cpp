#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

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

void RunExperiment(int m_k_cluster_mass[], const int size, const double treshold, const int seed)
{
    int* lattice= (int*)std::calloc(size*size, sizeof(int));
    PopulateLattice(lattice, size, treshold, seed);
    int k_cluster_label= 2;

    for (long int ij= 0; ij < size*size; ij++)
        if (lattice[ij] == 1)
            CheckSite(lattice, m_k_cluster_mass, k_cluster_label, size, ij);        
    std::free(lattice);
}

int main(int argc, char* argv[])
{
    if (argc < 5) // check command line arguments and give some help
     {  
        std::cerr << "Usage: " << argv[0]
                  << " int_occupation_prob_steps int_lattice_side_length int_initial_seed int_num_samples \n" << std::endl;
        return 1;
     }
    const int treshold_steps= atoi(argv[1]);
    const int size= atoi(argv[2]);
    const int initial_seed= atoi(argv[3]);
    const int num_samples= atoi(argv[4]);

    const double treshold_step_size= 1.0 / treshold_steps; // Build distribution of occupation probs
    std::vector<double> tresholds{0.592746}; // Init vector of probs and make sure p_c is included
    for (auto i= 1; i <= treshold_steps; i++)
        tresholds.push_back(i * treshold_step_size);

    for (auto& treshold : tresholds) // Run experiments on all chosen occupation probs
    {
        double histogram[size*size + 1]{0.0};
        std::cout << "#p=" << treshold << std::endl;
        for (auto n= 0; n < num_samples; n++) // Run num_samples amount of experiments
        {
            int* m_k_cluster_mass= (int*)std::calloc(size*size, sizeof(int));
            RunExperiment(m_k_cluster_mass, size, treshold, initial_seed + n);
            for (int k= 2; k < size*size; k++) // Record this experiment run's result in histogram
                if (m_k_cluster_mass[k] > 0)
                    histogram[m_k_cluster_mass[k]] += 1.0; 
            std::free(m_k_cluster_mass);
        }
        for (auto h= 1; h <= size*size; h++)
        {
            histogram[h] /= (num_samples * size * size); // Normalize histogram
            if (histogram[h] > 0)
                std::cout << h << "\t" << histogram[h] << std::endl; // Print if bin is non-zero
        }
        std::cout << "\n\n";
    }
    return 0;
}