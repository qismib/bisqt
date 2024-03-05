import math
import numpy as np

from qiskit.quantum_info import Pauli, Operator
from qiskit.quantum_info import SparsePauliOp

#Bose-Hubbard Hamiltonian of a system of bosons 
#binary mapping

class BoseHubbardHamiltonian: 
        
    def __init__(self, n_sites, single_site_qubits, geometry):
        
        # n_sites (int) = number of sites in the lattice (n_sites >= 2)
        # single_site_qubits (int) = number of qubits used to describe each site
        # geometry (LineLattice or another type of lattice) = how the sites are distributed in space and their connectivity 
              
        self.n_sites = n_sites
        self.ss_q = single_site_qubits
        self.geometry = geometry.to_adjacency_matrix() # the adjacency matrix gives informations about the connectivity of the sites

        # call methods to initialize attributes
        self.second_quantization_ops()
        self.lattice_connectivity()
        self.kinetic_energy_op()
        self.potential_energy_op()

    
    #creation, annihilation, number and identity operators on single_site_qubits number of qubits
    def second_quantization_ops(self):
       
        dim =  2**(self.ss_q)
        create_matrix = np.zeros((dim,dim),dtype=float)  
        
        for i in range(dim-1): 
            create_matrix[i+1][i] = math.sqrt(i+1)
        
        I_q_matrix = np.eye(dim, dtype=float)
        I_q_Op = Operator(I_q_matrix)
        
        # from matrix to Operator and from Operator to SparsePauliOp
        createOp = Operator(create_matrix)
        self.create = SparsePauliOp.from_operator(createOp)
     
        self.annihilate = self.create.adjoint()
        self.number = self.create.compose(self.annihilate).simplify()

        self.I_q =  SparsePauliOp.from_operator(I_q_Op)
        
    #removal of the diagonal terms in the adjacency matrix of the lattice
    def lattice_connectivity(self):

        for i in range(self.n_sites):
            for j in range(self.n_sites):
                if i == j:
                    self.geometry[i][j] = 0
        

    #kinetic energy operator
    def kinetic_energy_op(self): 

        Op = []
        for i in range(self.n_sites): 
            Op.append(self.I_q)

        kin = []
        for i in range(self.n_sites):
            for j in range(self.n_sites):
       
                if self.geometry[i][j] == 1:
                    Op[i] = self.create
                    Op[j] = self.annihilate
                      
                    Operator = Op[0]
                    for k in range (1 , n_sites ) :
                        Operator = Operator^Op[k]

                    kin.append(Operator)

                    Op[i] = self.I_q
                    Op[j] = self.I_q
            

        self.Kinetic = sum(kin).simplify()
        return self.Kinetic
              
    #potential energy (total interaction term) operator    
    def potential_energy_op(self):
        
        I_n_list = []
        Op = []
        for i in range(self.n_sites):
            I_n_list.append(self.I_q)
            Op.append(self.I_q)

        #identità per n siti con ss_q qubit per sito
        I_n = I_n_list[0]
        
        for k in range (1 , n_sites ) :
            I_n = I_n^I_n_list[i]
            
        
        self.I_n = I_n
        
        number_tot = []
        potent = []
        
        for i in range(self.n_sites): 
            Op[i] = self.number
       
            Operator = Op[0]
            
            for k in range (1 , n_sites ) :
                Operator = Operator^Op[k]

            number_tot.append(Operator)  #list of single site number operators on n_sites*ss_q qubit
                    
            Op[i] = self.I_q

        self.total_number = sum(number_tot)

        for i in range(self.n_sites):
            potent.append( number_tot[i].compose(number_tot[i] - self.I_n))


        self.Potential = sum(potent).simplify()

        return self.Potential, self.total_number
    
    #additional term in the Hamiltonian that allows to seletect states with a miximum number of particles
    # n = possible number of particles 
    def eigenvalues_selection_op(self, A, N):
        self.constraint_op = A*((N * self.I_n - self.total_number).compose((N * self.I_n - self.total_number)))     
        return self.constraint_op
    
    def get_H(self, t, U):
             
        self.U = U
        self.t = t
        H =  t*self.Kinetic + U*self.Potential
        return H
        