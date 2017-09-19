import numpy as np


class espack(object):
    def __init__(self, x, y, z, r):
        self.x = x  #length of carton domain
        self.y = y  #Width of carton domain
        self.z = z  #Height of carton domain
        self.r = r  #fruit radius
        
    # Used by stagger to calculate the number of fruit in a box height (for staggered)
    def closest_y_gap(self):
        self.node_number_y_temporary = np.ceil((self.y - (2 * self.r))/(self.r))     #initiallizes the maximum number of node values along y-axis
        for i in np.arange(self.node_number_y_temporary + 1):     #for each the possible nodes across y-axis  
            self.node_gap_y = ( self.y - ( 2 * self.r ) ) / (self.node_number_y_temporary - i - 1)     #calculates gap space (starts small) and itterates to bigger
            if np.sqrt(self.node_gap_x**2 + self.node_gap_y**2) > self.r*2:     #test to see if fruit are in contact
                self.node_number_y = self.node_number_y_temporary - i     #If not, it returns the number of spaces        
                break
                
        """
        node_number_x - Gaps along the x axis (width)
        node_gap_x - Space between x gaps (width)
        node_number_y - Gaps along the y axis (length)
        node_gap_y - Space between y gaps (length)
        """       
        
    def stagger(self):
        self.grid = np.zeros([4])     #initialise array for grid
        
        for i_x in np.arange(2,   (np.floor((self.x - (2 * self.r)) / (self.r)))   + 1):     #iterate through the number of possible grid nodes in the x-axis
            self.node_number_x = i_x
            self.node_gap_x = ( self.x - ( 2 * self.r ) ) / (i_x - 1)
            self.closest_y_gap()     #Find the most number of y-axis nodes without fruit contact   
            
            for i_y in np.arange(2, self.node_number_y + 1):     #iterate through less and less nodes along y-axis
                self.node_gap_y = ( self.y - ( 2 * self.r ) ) / (i_y - 1)
                self.node_gap_z = max(np.sqrt( abs(((2 * self.r)**2) - self.node_gap_x**2) ), np.sqrt( abs(((2 * self.r)**2) - self.node_gap_y**2) ) )
                
                if self.node_gap_z <= self.r:
                    self.node_gap_z = self.r     #just make it r if trays are too close together
                self.node_number_z = np.floor((self.z - 2 * self.r) / self.node_gap_z)                          #Find number of nodes across z-axis
                
                self.point = np.array([ i_x, i_y, self.node_number_z, i_x * i_y * self.node_number_z ])         #Generate nodes for every possible/feasable grid of nodes
                self.grid = np.vstack((self.grid, self.point))                                                  #Compile full list grid / (grid[np.argmax(grid, axis=0)[3]]) = select optimal grid array

        self.opt_number = self.grid[np.argmax(self.grid, axis=0)[3]]                                             #selction of optimal grid
        self.opt_gap = np.array([  ((self.x - 2 * self.r)/self.opt_number[0]) , (((self.y - 2 * self.r)/self.opt_number[1])) , (self.node_gap_z)  ])  # respective gap sizes
        
        #create an array of grid
        self.fruit_pos_x = np.linspace(self.r, self.opt_number[0]*self.opt_gap[0] + self.r, self.opt_number[0])  #Positions / no overlap
        self.fruit_pos_y = np.linspace(self.r, self.opt_number[1]*self.opt_gap[1] + self.r, self.opt_number[1])  #Positions / no overlap
        self.fruit_pos_z = np.linspace(self.r, self.z - self.r, self.opt_number[2])                              #Positions / fruit do overlap
                
        # Determine coordinates using gaps between nodes
        self.coordinates  = np.zeros([3])
        self.fruit_no = 0
        for i_z in np.arange(len(self.fruit_pos_z)):
            for i_x in np.arange(len(self.fruit_pos_x)):
                for i_y in np.arange(abs((i_z % 2) - (i_x % 2)), len(self.fruit_pos_y), 2):
                    
                    self.fruit_no = self.fruit_no + 1
                    
                    self.point = np.array([ self.fruit_pos_x[np.int(i_x)], self.fruit_pos_y[np.int(i_y)], self.fruit_pos_z[np.int(i_z)] ])
                    self.coordinates  = np.vstack((self.coordinates, self.point))
        
        # Coordinates of each fruit, without the initial (0,0,0) node
        self.coordinates  = self.coordinates[1:]
        
        # Volume of fruit in box
        self.fruit_vol = self.fruit_no*((4/3) * np.pi * self.r**3)
        
        # Volume of empty box
        self.box_vol = self.x * self.y * self.z
        
        # Void fraction of packed carton (porosity)
        self.por = (self.box_vol - self.fruit_vol) / self.box_vol
        
        # Fruit density in carton (porosity)
        self.den = 1 - self.por
        
        
    # simple cubic lattice
    def SC(self):
        self.node_number_x = np.int(self.x/(self.r*2) )
        self.node_number_y = np.int(self.y/(self.r*2) )
        self.node_number_z = np.int(self.z/(self.r*2) )
        
        self.coordinates = np.zeros([3])
        for i in range (self.node_number_x + 1):
            for j in range (self.node_number_y + 1):
                for k in range (self.node_number_z + 1):
                    self.el_x = i*self.r*2 + self.r
                    self.el_y = j*self.r*2 + self.r
                    self.el_z = k*self.r*2 + self.r
                    if self.el_x < (self.x - self.r) and self.el_y < (self.y - self.r) and self.el_z < (self.z - self.r) and self.el_x >= 0 and self.el_y >= 0 and self.el_z >= 0:
                    
                        self.point = np.array([self.el_x, self.el_y, self.el_z])
                        self.coordinates = np.vstack((self.coordinates, self.point))
                    
        # Coordinates of each fruit, without the initial (0,0,0) node                
        self.coordinates = self.coordinates[1:]
        
        # Volume of fruit in box
        self.fruit_vol = self.coordinates.shape[0]*((4/3) * np.pi * self.r**3)
        
        # Volume of empty box
        self.box_vol = self.x * self.y * self.z
        
        # Void fraction of packed carton (porosity)
        self.por = (self.box_vol - self.fruit_vol) / self.box_vol
        
        # Fruit density in carton (porosity)
        self.den = 1 - self.por 
        
    
    # Hexagonal close packing python
    def HPC(self):
        
        self.node_number_x = np.int(self.x/self.r)
        self.node_number_y = np.int(self.y/self.r)
        self.node_number_z = np.int(self.z/self.r)
        
        self.coordinates = np.zeros([3])

        for i in range(self.node_number_x + 1):
            for j in range(self.node_number_y + 1):
                for k in range(self.node_number_z + 1):
                    self.el_x = (2*i + ((j + k)%2))*self.r + self.r
                    self.el_y = (np.sqrt(3)*(j + (1/3)*(k%2)))*self.r + self.r
                    self.el_z = (((2*np.sqrt(6))/3)*k)*self.r + self.r
                    
                    if self.el_x < (self.x - self.r) and self.el_y < (self.y - self.r) and self.el_z < (self.z - self.r) and self.el_x >= 0 and self.el_y >= 0 and self.el_z >= 0:
                        self.point = np.array([self.el_x, self.el_y, self.el_z])
                        self.coordinates = np.vstack((self.coordinates, self.point))
                        
        # Coordinates of each fruit, without the initial (0,0,0) node                
        self.coordinates = self.coordinates[1:]
        
        # Volume of fruit in box
        self.fruit_vol = self.coordinates.shape[0]*((4/3) * np.pi * self.r**3)
        
        # Volume of empty box
        self.box_vol = self.x * self.y * self.z
        
        # Void fraction of packed carton (porosity)
        self.por = (self.box_vol - self.fruit_vol) / self.box_vol
        
        # Fruit density in carton (porosity)
        self.den = 1 - self.por
    
    
    # Body-centered cubic (BCC)
    def BCC(self):
        self.coordinates = np.zeros([3])
        self.root = (np.sqrt(3/4))
        self.maxx = np.int(np.amax ([self.x/self.r,self.y/self.r,self.z/self.r]))
        
        for i in range(self.maxx):
            for j in range(self.maxx):
                for k in range(self.maxx):
                    self.el_x = 2*self.r*((0.5*(-i+j+k))/self.root) + self.r
                    self.el_y = 2*self.r*((0.5*(i-j+k))/self.root) + self.r
                    self.el_z = 2*self.r*((0.5*(i+j-k))/self.root) + self.r
                    if self.el_x < (self.x-self.r) and self.el_y < (self.y-self.r) and self.el_z < (self.z-self.r) and self.el_x >= 0 and self.el_y >= 0 and self.el_z >= 0:
                        self.point = np.array([self.el_x, self.el_y, self.el_z])
                        self.coordinates = np.vstack((self.coordinates, self.point))
                        
        # Coordinates of each fruit, without the initial (0,0,0) node  
        self.coordinates = self.coordinates[1:]
        
        # Volume of fruit in box
        self.fruit_vol = self.coordinates.shape[0]*((4/3) * np.pi * self.r**3)
        
        # Volume of empty box
        self.box_vol = self.x * self.y * self.z
        
        # Void fraction of packed carton (porosity)
        self.por = (self.box_vol - self.fruit_vol) / self.box_vol
        
        # Fruit density in carton (porosity)
        self.den = 1 - self.por         
        
        
    # Face-centered cubic close packing
    def FCC(self):
        self.coordinates = np.zeros([3])
        self.root = np.sqrt(0.5)
        self.maxx = np.int(np.amax ([self.x/self.r,self.y/self.r,self.z/self.r]))
        for i in range(-self.maxx, self.maxx):
            for j in range(-self.maxx, self.maxx):
                for k in range(-self.maxx, self.maxx):
                    self.el_x = 2*self.r* ((0.5*(j+k))/self.root) + self.r
                    self.el_y = 2*self.r* ((0.5*(i+k))/self.root) + self.r
                    self.el_z = 2*self.r* ((0.5*(i+j))/self.root) + self.r
                    if self.el_x < (self.x-self.r) and self.el_y < (self.y-self.r) and self.el_z < (self.z-self.r) and self.el_x >= 0 and self.el_y >= 0 and self.el_z >= 0:
                        self.point = np.array([self.el_x, self.el_y, self.el_z])
                        self.coordinates = np.vstack((self.coordinates, self.point))
                        
        # Coordinates of each fruit, without the initial (0,0,0) node  
        self.coordinates = self.coordinates[1:]
        
        # Volume of fruit in box
        self.fruit_vol = self.coordinates.shape[0]*((4/3) * np.pi * self.r**3)
        
        # Volume of empty box
        self.box_vol = self.x * self.y * self.z
        
        # Void fraction of packed carton (porosity)
        self.por = (self.box_vol - self.fruit_vol) / self.box_vol
        
        # Fruit density in carton (porosity)
        self.den = 1 - self.por  
        
        
