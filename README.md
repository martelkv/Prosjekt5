# Disease Modeling    
*Authors: Jelena Jankovic and Marte Lunde Kvam*. 

### How to run the program
- ```g++ *.cpp -o a``` to compile.  
- ```./a``` to run.  

### How to obtain the plots
To obtain the plots, first run the program, then execute  ```python plot.py```.  
  - This command will run the file **plot.py** and use the txt-files in the **output/**-folder.  
To obtain the plots from the results-section in the report, change these lines to *True* in **plot.py** like this:  
```python3
# To see the plots from section 3.1
a = True

# To see the plots from section 3.2
b = True

# To see the plots from section 3.3
c_MonteCarlo = True

# To see the plots from section 3.4
d_MonteCarlo = True

# To see the plots from section 3.5
e_MonteCarlo = True

# To obtain the results from the standard deviation- and expectation values
std = True
std_rungakutta= True
mean = True

```

### Requirements
To be able to run the programs, you need the following packages installed:  
- Numpy    
- Armadillo  
- c++  
- matplotlib   
