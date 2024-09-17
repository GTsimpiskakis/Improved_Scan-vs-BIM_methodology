#Point cloud and IFC overlaping (matching)
#This python script takes as an input the resulted segmented .ply file which is sufficiently downsampled 
#and now can be analysed in Blender and transformed efficiently
# The user should manually capture the coordinates of the points A1 in IFC file, and O2 and A2 in ,ply file.
#as these coordinates will be used as inputs in below transformation process.

import open3d as o3d
import numpy as np
import math
import os 
# Construct the rotation matrix for rotation around x-axis
def rotation_about_x_axis(theta):
    rotation_matrix = np.array([
    [1, 0, 0, 0],
    [0, np.cos(theta), np.sin(theta), 0], 
    [0, -np.sin(theta), np.cos(theta), 0],
    [0, 0, 0, 1]
    ])
    #array([[1, 0, 0, 0.],
           #[0, cos(θ), sin(θ), 0.],
           #[0, -sin(θ), cos(θ), 0.],
           #[0., 0., 0., 1.]])
    return rotation_matrix
# Construct the rotation matrix for rotation around y-axis
def rotation_about_y_axis(theta):
    rotation_matrix = np.array([
    [np.cos(theta), 0, np.sin(theta), 0],
    [0, 1, 0 , 0], 
    [np.sin(theta), 0, np.cos(theta), 0],
    [0, 0, 0, 1]
    ])
    #array([[cos(θ), 0, sin(θ), 0.],
           #[0, 1, 0., 0.],
           #[sin(θ), 0., cos(θ)., 0.],
           #[0., 0., 0., 1.]])
    return rotation_matrix

def rotation_about_z_axis(theta):
    # Construct the rotation matrix for rotation around z-axis
    rotation_matrix = np.array([
    [np.cos(theta), -np.sin(theta), 0, 0],
    [np.sin(theta), np.cos(theta), 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
    ])
    #array([[cos(θ), -sin(θ), 0., 0.],
           #[sin(θ),cos(θ), 0., 0.],
           #[0., 0., 1., 0.],
           #[0., 0., 0., 1.]])
    return rotation_matrix
# Create a point cloud
pcd_file_path = r"C:/Python Output/Real Data/PLY files/SDK_real_dataset_final/PLY_SDK_EU_Real_Merged_PCM_Translated.ply"
directory, file_name_with_extension = os.path.split(pcd_file_path)
file_name, extension = os.path.splitext(file_name_with_extension)
pcd = o3d.io.read_point_cloud(pcd_file_path)
translation_matrix = np.identity(4) #this command creates the translation matrix M, se below
    #array([[1., 0., 0., 0.],
    #[0., 1., 0., 0.],
    #[0., 0., 1., 0.],
    #[0., 0., 0., 1.]])
Transformation_type =int(input("What operation do you want to do? For translation press 1, For Rotation press 2 for both press 3"))
if Transformation_type == 1: #Then Only translation
    print("Give the absolute Coordinates in meters (m) of the O2 (point cloud origin point) in terms of O1(0,0,0) (IFC model origin point) ")
    x = float(input("Give cooridnate x in terms of the current origin point"))
    y = float(input("Give cooridnate y in terms of the current origin point"))
    z = float(input("Give cooridnate z in terms of the current origin point"))
    # Define translation vector
    translation = np.array([-x, -y, -z]) #this command specifies the translation vector M2=[tx,ty,tz,1]
    # Create transformation matrix for translation
    translation_matrix[:3, 3] = translation #this command adds the 3x1 translation vector M2 in the 
    #translation matrix M, therefore gives values to the tx, ty,tz.
    #array([[1., 0., 0., 1.],
    #[0., 1., 0., 0.],
    #[0., 0., 1., 0.],
    #[0., 0., 0., 1.]])
    # Define two vectors
    rotation_matrix = np.array([
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1]
    ])
    additional_word = "_Translated"
    # Modify the filename based on the condition
    new_file_name = file_name + additional_word
    # Recreate the output file paths
    pcd_file_path_out = [os.path.join(directory, new_file_name + extension)]
    #pcd_file_path_out = r"C:/Python Output/Testdata/PLY files/PLY_GT_TOPOPLAN_TEST1 1_Transl_Rotated Form A2 to A1_test 3.ply"
elif Transformation_type ==2 or Transformation_type ==3 : # only Rotation or Both Translation and Rotation
    if Transformation_type ==3:#Both Translation and Rotation simultaneously
        print("Give Coordinates of the preferable O(0,0,0) origin point in point cloud")
        x = float(input("Give cooridnate x in terms of the current origin point"))
        y = float(input("Give cooridnate y in terms of the current origin point"))
        z = float(input("Give cooridnate z in terms of the current origin point"))
        # Define translation vector
        translation = np.array([-x, -y, -z]) #this command specifies the translation vector M2=[tx,ty,tz,1]
        # Create transformation matrix for translation
        translation_matrix[:3, 3] = translation #this command adds the 3x1 translation vector M2 in the 
        #translation matrix M, therefore gives values to the tx, ty,tz.
        #array([[1., 0., 0., 1.],
        #[0., 1., 0., 0.],
        #[0., 0., 1., 0.],
        #[0., 0., 0., 1.]])
        # Define two vectors
        additional_word = "_Translated_Rotated"
        # Modify the filename based on the condition
        new_file_name = file_name + additional_word
        # Recreate the output file paths
        pcd_file_path_out = [os.path.join(directory, new_file_name + extension)]
       
    else:
        additional_word = "_Rotated1"
        # Modify the filename based on the condition
        new_file_name = file_name + additional_word
        # Recreate the output file paths
        pcd_file_path_out = [os.path.join(directory, new_file_name + extension)]
    #Spcify the points     
    print("Give Coordinates of the image Point A1") # For this example: A1(0.5,0,0)
    x1 = float(input("Give cooridnate x in terms of the current origin point"))
    y1 = float(input("Give cooridnate y in terms of the current origin point"))
    z1 = float(input("Give cooridnate z in terms of the current origin point"))
    print("Give Coordinates of the Point A2, which should be transformed to A1")#For this example: A2(0.433,0.258,0)
    x2 = float(input("Give cooridnate x in terms of the current origin point"))
    y2 = float(input("Give cooridnate y in terms of the current origin point"))
    z2 = float(input("Give cooridnate z in terms of the current origin point"))

    #set the Z coordinates to zero in order to convert the vectors in XY plane
    if z1 !=0: 
        z1 = 0
    if z2 !=0:
        z2 = 0
    v_A1 = np.array([x1, y1, z1])
    v_A2 = np.array([x2, y2, z2])

    # Calculate the dot product
    dot_product = np.dot(v_A2, v_A1) #We calculate the angle from A2 to A1, Be careful this is different from the previous examples in Code 12,
    #in previous examples we had:dot_product = np.dot(v_A1, v_A2), namely angle from A1 to A2.
    cross_product = np.cross(v_A2, v_A1)
    Cross_Product_sign_z3 = cross_product[2]
    
    # Calculate the magnitude
    magn_A1 = np.linalg.norm(v_A1)
    magn_A2 = np.linalg.norm(v_A2)

    print("Magnitudes A1,A2:", magn_A1, magn_A2)
    print("Dot product:", dot_product)
    print("Cross product:", cross_product)
    print("Cross_Product_sign_z value:", Cross_Product_sign_z3)
    # Cosine value
    cos_theta = dot_product/(magn_A1*magn_A2)

    # Calculate the arccos of the cosine value, this theta is the counterclock wise theta
    theta = math.acos(cos_theta)
    print("Theta (in radians):", theta)
    # Convert radians to degrees only for prinitng purpose
    theta_degrees = math.degrees(theta)
    print("Theta (in degrees):", theta_degrees)
        #as the cos(θ)=Α1.Α2/|Α1|x|A2| (Equation 1) always produces the angle between A1 and A2 that is lower than 180 degrees.
        #We always want the angle that goes clockwise from A2 to A1 which deemed negative in python environment.
        #thus, we investigate in Blender software and we see if the desired angle is below or above 180, this can easily been seen by a naked eye.
    #if the angle is above 180 degrees then the equation will produce the counterclockwise angle from A2 to A1, but we always want the clockwise angle 
    #from A2 to A1. That is why we substracted from 360 degrees only if the deisre angle is above 180 degrees.
    if Cross_Product_sign_z3 >= 0:
        theta_clockwise = theta # this means that the deisre angle is below 180 degrees so it is produced correctly by the equation 1 cos(θ)
        print(theta_clockwise)
        
        #we give a minus sign in θ' (-θ') to illustrate the clockwise rotation (See diagrams 06/03/2024)
    else:
        theta_clockwise = 2*(np.pi) - theta # θ'=360-θ
        print(theta_clockwise)
    r= int(input("give the axis of rotation, 1 for x, 2 for y, 3 for z"))
    #Specify the rotation matrix
    if r==1:
        rotation_matrix = rotation_about_x_axis(theta_clockwise)
    elif r ==2:
        rotation_matrix = rotation_about_y_axis(theta_clockwise)
    elif r ==3:
        rotation_matrix = rotation_about_z_axis(theta_clockwise)
# Combine translation and rotation matrices
combined_matrix = np.dot(translation_matrix, rotation_matrix)

# Apply translation to the point cloud
pcd.transform(combined_matrix)
def write_ply_from_pcd(points, colors, filename):
    vertex = np.zeros(points.shape[0], dtype=[("x", "f4"), ("y", "f4"), ("z", "f4"), ("red", "u1"), ("green", "u1"), ("blue", "u1")])
    vertex["x"] = points[:, 0]
    vertex["y"] = points[:, 1]
    vertex["z"] = points[:, 2]
    vertex["red"] = colors[:, 0]
    vertex["green"] = colors[:, 1]
    vertex["blue"] = colors[:, 2]

    with open(filename, 'w') as f:
        f.write("ply\n")
        f.write("format ascii 1.0\n")
        f.write("element vertex {}\n".format(len(points)))
        f.write("property float x\n")
        f.write("property float y\n")
        f.write("property float z\n")
        f.write("property uchar red\n")
        f.write("property uchar green\n")
        f.write("property uchar blue\n")
        f.write("end_header\n")
        for i in range(len(points)):
            f.write("{:.3f} {:.3f} {:.3f} {:.0f} {:.0f} {:.0f}\n".format(points[i][0], points[i][1], points[i][2], colors[i][0], colors[i][1], colors[i][2]))
points = np.asarray(pcd.points)
colors = np.asarray(pcd.colors)
write_ply_from_pcd(points, colors*255, pcd_file_path_out[0] )
# Visualize original and translated point clouds
o3d.visualization.draw_geometries([pcd])

