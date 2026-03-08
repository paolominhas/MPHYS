import uproot
import json

# 1. Open the Geant4 root file
file = uproot.open("../../initial-proj/final_sim_data/250proton.root")

# 2. Access the 'Tree' (the table) that contains your particle data
# Change "Events" to whatever your TTree is named
tree = file["hibeam"] 

# 3. Extract the arrays (Change these string names to match your Geant4 branches)
# We are assuming you saved the origin (x,y,z) and momentum/velocity (px,py,pz)
x_data = tree["target_Position_X"].array(library="np")
y_data = tree["target_Position_Y"].array(library="np")
z_data = tree["target_Position_Z"].array(library="np")
px_data = tree["target_Momentum_X"].array(library="np")
py_data = tree["target_Momentum_Y"].array(library="np")
pz_data = tree["target_Momentum_Z"].array(library="np")
pdg_codes = tree["target_PDG"].array(library="np") # For color coding

# 4. Format it into a clean JSON structure
# 4. Format it into a clean JSON structure
particles = []
global_particle_id = 0 # Keeps track of the ID across all events

# Outer loop: Go through every event
for i in range(len(x_data)):
    
    # Inner loop: Go through every particle inside the current event
    for j in range(len(x_data[i])):
        
        # Extract the specific PDG code for this one particle
        pdg = int(pdg_codes[i][j])
        
        # Assign colors based on particle type
        color = "#ffffff" # Default white
        if pdg == 11 or pdg == -11: 
            color = "#00ffff" # Cyan (Electrons)
        elif pdg == 13 or pdg == -13: 
            color = "#ff00ff" # Magenta (Muons)
        elif pdg == 2212: 
            color = "#ff3333" # Red (Protons - seen in your data!)
        elif pdg > 1000000000: 
            color = "#ffff00" # Yellow (Heavy Nuclei / Isotopes - seen in your data!)

        # Extract the exact coordinates and momentum for this particle
        particle = {
            "id": global_particle_id,
            "type": pdg,
            "origin": [float(x_data[i][j]), float(y_data[i][j]), float(z_data[i][j])],
            "velocity": [float(px_data[i][j]), float(py_data[i][j]), float(pz_data[i][j])],
            "color": color
        }
        
        particles.append(particle)
        global_particle_id += 1

output_data = {
    "dataset_name": "250proton",
    "total_particles": len(particles),
    "particles": particles
}

# 5. Save it to a JSON file
with open("collision-data.json", "w") as outfile:
    json.dump(output_data, outfile, indent=2)

print("YAY successfully extracted ROOT data to collision-data.json!")