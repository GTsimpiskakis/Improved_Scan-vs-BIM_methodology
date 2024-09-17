import ifcopenshell
import trimesh
import open3d as o3d
import numpy as np
import math
import os 
import re
import pandas as pd
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas

#ifc_file_path = r"C:/Python Output/Testdata/IFC files/GT_Topoplan_psich_column1_More columns_Pipes.ifc"
ifc_file_path = r"C:/Python Output/Real Data/IFC files/Eurostate_Tower_04_floor_Final_060624.ifc"
# Open the IFC file
ifc_file = ifcopenshell.open(ifc_file_path)
# Get all desirable objects
slabs = ifc_file.by_type("IfcSlab")
Columns = ifc_file.by_type("IfcColumn")
walls = ifc_file.by_type("IfcWall")
beams = ifc_file.by_type("IfcBeam")
pipes = ifc_file.by_type("IfcPipeSegment")
# Combine all objects into a single numpy array
IfcObjects= walls + Columns + beams + slabs + pipes
O = 0., 0., 0.
X = 1., 0., 0.

#this is a vector that defines the positive direction of x axis
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
translation_matrix = np.identity(4) #this command creates the translation matrix M, se below
    #array([[1., 0., 0., 0.],
    #[0., 1., 0., 0.],
    #[0., 0., 1., 0.],
    #[0., 0., 0., 1.]])
def object_ObjectPlacement_inheritence(object):
  
    print(object)
    objectPlacement1 = object.ObjectPlacement
    print(f"LocalPlacement1",objectPlacement1) #First IfcLocalPlacement which refers to wall
    PlacementRelTo1 = objectPlacement1.PlacementRelTo #second IfcLocalPlacement which referes to the IfcBuildingstorey
    print(f"PlacementRelTo1",PlacementRelTo1) #second IfcLocalPlacement which referes to the IfcBuildingstorey
    RelativePlacement1 = objectPlacement1.RelativePlacement #this produces the IfcAxis2Placement3D which referes to the wall 
    print(f"RelativePlacement1:",RelativePlacement1)
    #For first IfcAxis2Placement3D
    print("for first IfcAxis2Placement3D")
    Location1 = RelativePlacement1.Location #this produces the Transforamtion needed for the wall origin point in 
    #order to match with the Building storeys origin point
    print(f"Location1:",Location1)
    Local_Origin_Ifcobject = Location1.Coordinates
    print(f"Local_Origin_Ifcobject Coordinates:",Local_Origin_Ifcobject)
    Axis1 = RelativePlacement1.Axis
    print(f"Axis:",Axis1)
    RefDirection1 = RelativePlacement1.RefDirection
    print(f"RefDirection1:",RefDirection1)
    if RefDirection1 is not None:
        #Then the SweptArea of this IfcElement is IfcArbitraryClosedProfileDef namely its LCS is rotated in terms of GCS
        DirectionRatios = RefDirection1.DirectionRatios
        print(f"DirectionRatios:",DirectionRatios)
    else:
        DirectionRatios = None
    #In this part we analyse the second IfcLocalPlacement which refers to IfcBuildingStorey.
    print("second ifcLocalPlacement")
    print(f"Analyse PlacementRelTo1",PlacementRelTo1) 
    PlacementRelTo2= PlacementRelTo1.PlacementRelTo
    print(f"PlacementRelTo2",PlacementRelTo2) #produces a third IfcLocal Placement which refers to IfcBuilding
    RelativePlacement2= PlacementRelTo1.RelativePlacement
    print(f"RelativePlacement2",RelativePlacement2) #Produces an IfcAxis2Placement3D which refers to Ifc BuildingStorey
    #For second IfcAxis2Placement3D
    Location2 = RelativePlacement2.Location #From IfcAxis2Placement3D we find the necessary transformation needed 
    #for an IfcBuildingStorey in order its coordinations match with the CRS of its parent's namely IfcBuilding.
    print(f"Location2:",Location2)

    #In this part we analyse the third IfcLocalPlacement which refers to IfcBuilding
    print("third ifcLocalPlacement")
    print(f"Analyse PlacementRelTo2",PlacementRelTo2)
    PlacementRelTo3= PlacementRelTo2.PlacementRelTo
    print(f"PlacementRelTo3",PlacementRelTo3) #produces a Fourth IfcLocal Placement which refers to IfcSite
    RelativePlacement3 = PlacementRelTo2.RelativePlacement
    print(f"RelativePlacement2:",RelativePlacement3) #Produces an IfcAxis2Placement3D which refers to IfcBuilding
    #For second IfcAxis2Placement3D
    Location3 = RelativePlacement3.Location #From IfcAxis2Placement3D we find the necessary transformation needed 
    #for an IfcBuilding in order its coordinations match with the CRS of its parent's namely IfcSite.
    print(f"Location3:",Location3)

    #In this part we analyse the Fourth IfcLocalPlacement which refers to IfcSite
    print("Fourth local placement")
    print(f"Analyse PlacementRelTo3",PlacementRelTo3)
    PlacementRelTo4 = PlacementRelTo3.PlacementRelTo 
    print("PlacementRelTo4", PlacementRelTo4)#produces a None value because IfcSite does not inherit its CRS
    #from any other parent entity
    RelativePlacement4 = PlacementRelTo3.RelativePlacement #Produces an IfcAxis2Placement3D which refers to IfcSite
    print(f"RelativePlacement3:",RelativePlacement4)
    #For second IfcAxis2Placement3D
    Location4 = RelativePlacement3.Location #Produces the origin point in which the IfcSite CRS it begins
    print(f"Location4:",Location4)

    return Local_Origin_Ifcobject, RefDirection1, DirectionRatios, Axis1, Location2

def object_Representations(object):
            #Identify object's "Representations", namely the IfcShapeRepresentations entities
            Representation = object.Representation
            print(Representation)
            Representations = Representation.Representations
            print(Representations)
            return Representations
def IfcGeometricRepresentationContext(Representations):
            ContextOfItems1 = Representations[0].ContextOfItems
            print(f"ContextOfItems1", ContextOfItems1)
            ParentContext1 = ContextOfItems1.ParentContext
            print(f"ParentContext1", ParentContext1)
            WorldCoordinateSystem = ParentContext1.WorldCoordinateSystem
            print(f"WorldCoordinateSystem", WorldCoordinateSystem)
            TrueNorth = ParentContext1.TrueNorth
            print(f"TrueNorth", TrueNorth)

def ShapeRepresentations_loop(Representations):
    #The Item IfcExtrudedAreaSolid which we are looking for is located in different locations inside the Ifc Documentation depending on the IfcElement
    #we are analyse. #For example: 
    #For IfcColumns the IfcExtrudedAreaSolid item is included inside the MappedItem and is found with the relationships: MappingSource-->MappedRepresentation-->Items"
    #For IfcWalls the IfcExtrudedAreaSolid item is included in the second ShapeRepresentation of the ProductDefinitionShpae
    #and is found with the relationship "Items" directly. That's why we use this for loop to iterate through any ShapeRepresentation 
    Item_IfcExtrudedAreaSolid = None
    Item_IfcIfcPolyline = None
    PolylinePoints = None
    for ShapeRepresentation in Representations:
        Item1 = ShapeRepresentation.Items
        print(Item1)
        if   re.search(patternMappedItem, str(Item1)):
            print("MappedItem Started")
            MappingTarget = Item1[0].MappingTarget
            print(f"MappingTarget",MappingTarget)
            MappingSource = Item1[0].MappingSource
            print(f"MappingSource",MappingSource)
            MappingOrigin = MappingSource.MappingOrigin
            print(f"MappingOrigin", MappingOrigin)
            MappedRepresentation = MappingSource.MappedRepresentation
            print(f"MappedRepresentation", MappedRepresentation)
            if re.search(patternIfcExtrudedAreaSolid , str(MappedRepresentation.Items)):
                Item_IfcExtrudedAreaSolid  = MappedRepresentation.Items
                print("The IfcExtrudedAreaSolid was found: ", Item_IfcExtrudedAreaSolid)
                ContextOfItems2 = MappedRepresentation.ContextOfItems
                print(f"ContextOfItems2", ContextOfItems2)
                
        elif re.search(patternIfcExtrudedAreaSolid , str(Item1)):
            Item_IfcExtrudedAreaSolid = Item1
            print("The IfcExtrudedAreaSolid was found: ", Item_IfcExtrudedAreaSolid)
        elif re.search(patternIfcPolyline, str(Item1)):
            #then we have Slab and, the First shape's representation Item represents theIfcExtrudedAreaSolid
            Item_IfcIfcPolyline = Item1
            print("The IfcPolyline was found: ", Item_IfcIfcPolyline)
            print("Footprint: Curved 2D Represent Identifier")
            PolylinePoints = Item1[0].Points
            print(f"PolylinePoints", PolylinePoints)

        elif re.search(patternClippingResultItem , str(Item1)):
            FirstOperand_item = Item1[0].FirstOperand
            print(f"IfcBooleanClippingResult", FirstOperand_item)
            Item_IfcExtrudedAreaSolid = FirstOperand_item
            print(f"Item_IfcExtrudedAreaSolid", Item_IfcExtrudedAreaSolid) 

    return Item_IfcExtrudedAreaSolid, Item_IfcIfcPolyline, PolylinePoints

def Item_IfcExtrudedAreaSolid_analysis(Item_IfcExtrudedAreaSolid):
    print("IfcExtrudedAreaSolid Started")
    if re.search(patternArbitary , str(Item_IfcExtrudedAreaSolid[0])) or re.search(patternRectangle , str(Item_IfcExtrudedAreaSolid[0])):
        SweptArea = Item_IfcExtrudedAreaSolid[0] #or SweptArea = Item_IfcExtrudedAreaSolid.SweptArea
        print(f"SweptArea", SweptArea)
        Position = Item_IfcExtrudedAreaSolid.Position
        print(f"Position", Position)
        ExtrudedDirection = Item_IfcExtrudedAreaSolid.ExtrudedDirection
        print(f"ExtrudedDirection", ExtrudedDirection)
        Depth = Item_IfcExtrudedAreaSolid.Depth
        print(f"Depth", Depth)
    else:
        SweptArea = Item_IfcExtrudedAreaSolid[0].SweptArea
        print(f"SweptArea", SweptArea)
        Position = Item_IfcExtrudedAreaSolid[0].Position
        print(f"Position", Position)
        ExtrudedDirection = Item_IfcExtrudedAreaSolid[0].ExtrudedDirection
        print(f"ExtrudedDirection", ExtrudedDirection)
        Depth = Item_IfcExtrudedAreaSolid[0].Depth
        print(f"Depth", Depth)
    print("For Axis2LocalPlacement3D of IfcExtrudedAreaSolid ")
    Location_of_IfcExtrudedAreaSolid = Position.Location
    print(f"Location_of_IfcExtrudedAreaSolid", Location_of_IfcExtrudedAreaSolid)
    Axis = Position.Axis
    print(f"Axis", Axis)
    RefDirection = Position.RefDirection
    print(f"RefDirection", RefDirection)
    return SweptArea, Position, ExtrudedDirection, Depth, Location_of_IfcExtrudedAreaSolid, Axis, RefDirection
def Item_SweptArea_Arbitary(SweptArea):
    ProfDefPoints = None
    print("Entity is IfcArbitraryClosedProfileDef. Continue the process.")
    OuterCurve = SweptArea.OuterCurve
    print(f"OuterCurve", OuterCurve)
    ProfDefPoints = OuterCurve.Points
    print(f"ProfDefPoints(IfcPolyline points)", ProfDefPoints)
    return ProfDefPoints

def Item_SweptArea_Rectangle(SweptArea):
    print("Entity is IfcRectangleProfileDef. Continue the process.")
    PositionRec = SweptArea.Position
    print(f"Position of Rectangle profile Definition", PositionRec)

    XDim = SweptArea.XDim
    print(f"XDim", XDim)
    YDim = SweptArea.YDim
    print(f"YDim", YDim)
    return XDim, YDim
def Step2_3_4_Vertices_Transformation(Vert_Hyp_LCS, x, y , z, LGC_direction_Vect):
    print("Step 2 Started: Translate the vertices from Hypothetical scenario in Oc to O origin point")
    #Step 2 "TRANSLATION": After we have deifned the Matrix for Ifc Element Vertices and Local origin point we should translate them to the Global Origin point
    #        Thus we transform Matrix: Vert_Hyp_LCS into Vert_Hyp_GCS by substracting the coordinates of the Oc(x,y,z)= Oc(3061,1786,0)
    #Calculations: translation_matrix([[1., 0., 0., Dx],
                                        #[0., 1., 0., Dy],
                                        #[0., 0., 1., Dz],
                                        #[0., 0., 0., 1.]])*(Ax,Ay,Az,1)= (Ax + Dx, Ay + Dy, Az + Dz, 1)
    # Define translation vector
    translation = np.array([-x, -y, -z]) #this command specifies the translation vector M2=[tx,ty,tz,1]
    # Create transformation matrix for translation
    translation_matrix[:3, 3] = translation #this command adds the 3x1 translation vector M2 in the 
    #translation matrix M, therefore gives values to the tx, ty,tz.
    #array([[1., 0., 0., Dx.],
    #[0., 1., 0., Dy.],
    #[0., 0., 1., Dz.],
    #[0., 0., 0., 1.]])
    # The Transformatino equation has certain structure as we explained before, thus, we need to transform point by point. This illustrated below.
    i=0
    Vert_Hyp_GCS = np.empty((0, Vert_Hyp_LCS.shape[1]))
    for i in range(8):
        Transformation1 = np.dot(translation_matrix,Vert_Hyp_LCS[i] )
        Vert_Hyp_GCS = np.vstack((Vert_Hyp_GCS, Transformation1))
        i = i+1
    print("The new matrix demonstrates the translated IfcElement vertices from the Hypothetical scenario in LCS to GCS",Vert_Hyp_GCS )
    print("Step 3 Started: Rotate the vertices in terms of z about angle θ in order GCS match with LCS")
    #Step 3 "ROTATION": a) Now we should rotate the new matrix in terms of angle θ. This angle defines the direction of the LCS.
    # So the translated matrix from LCS to GCS origin point will be rotated and a new matrix will be produce that will contain the coordinates of the 
    #vertices in the correct direction of LCS but in the origin point of GCS.
    # b) the angle θ can be found by the two normalised vectors (0-1) that IFC file provides for Ifcelement
    #       Vector1: GCS_x(1,0,0) is the direction of the positive X-axis of GCS.
    #       Vector2: LGC_direction_Vect(-0.9119, -0.4103, 0) is the direction of positive X-axis of LCS.

    print(" Coordinates of the image Point A1 Positive X-axis of GCS")
    #Here we store the normalise vector [0-1] that defines the direction of positive X-axis of GCS.
    GCS_x = np.array(X)
    x1 = GCS_x[0]
    y1 = GCS_x[1]
    z1 = GCS_x[2]

    print(" Coordinates of the Point A2, which should be transformed to A1")
    #Here we store the normalised vector [0-1] that defines the direction of positive x-axis of IfcElement's LCS.
    Rotated_vector = np.array(LGC_direction_Vect)
    x2 = Rotated_vector[0]
    y2 = Rotated_vector[1]
    z2 = Rotated_vector[2]
    #Regarding the two Vectors A1,A2, we calculate the angle θ, namely the angle that LCS should cover to be matched with GCS.
    #However, in this case we Rotate the GCS around Z-axis about angle θ in order to transform the coordinates of vertices in the direction of LCS
    v_A1 = np.array([x1, y1, z1])
    v_A2 = np.array([x2, y2, z2])

    # Calculate the dot product
    dot_product = np.dot(v_A1, v_A2)
    cross_product = np.cross(v_A1, v_A2)
    # Calculate the magnitude
    magn_A1 = np.linalg.norm(v_A1)
    magn_A2 = np.linalg.norm(v_A2)

    print("Magnitudes A1,A2:", magn_A1, magn_A2)
    print("Dot product:", dot_product)
    print("Cross product:", cross_product)
    Cross_Product_sign_z3 = cross_product[2]
    print("Cross_Product_sign_z3:", Cross_Product_sign_z3)
    # Cosine value
    cos_theta = dot_product/(magn_A1*magn_A2)
    print(cos_theta)

    # Calculate the arccos of the cosine value, this theta is the counterclock wise theta
    theta = np.arccos(cos_theta) #This command is always producing positive angle θ. that is why we always try to identify if it is more or less than 180 degrees
    print("Theta (in radians):", theta)
    # Convert radians to degrees only for prinitng purpose
    theta_degrees = math.degrees(theta) 
    print("Theta (in degrees):", theta_degrees)
    # in this phase we will always desire to find the angle θ that runs from vecrtor A1 to vector A2, This angle cannot be specified by the dot product
    #As the equation cos(θ)=Α1.Α2/|Α1|*|A2| does not take into account the direction of angle θ (From A1 to A2). In fact this equation always producing
    #the θ between A1 and A2 that is less than 180 degrees. Hence, when our desired rotation angle from A1 to A2 is bigger than 180 degrees the equation
    # will produce wrong θ result.
    #To solve that problem we exploit the properties of cross product between A1 and A2. The Equation 2: A3 = Α1xΑ2 = |Α1|*|A2|sin(θ) takes into account the
    #direction of angle θ which in this case it goes from A1 to A2. 
    # So if the angle θ from A1 to A2 is below 180 degrees then A3(x3,y3,z3)>0 namely, z3>0,
    #    if the angle θ from A1 to A2 is above 180 degrees then A3(x3,y3,z3)>00 namely, z3<0
    #So by calculating the cross product A3 = A1xA2 we can specify the correct θ (from Α1 to Α2) that is produced by the Equation 1 of dot product
    if Cross_Product_sign_z3>=0: #then θ from A1 to A2 is below 180 degrees so the Equation 1 of Dot product have produced it correctly
        theta_New = theta #Theta Counterclockwise
        
    else: #Cross_Product_sign_z3<0 which means that the θ from A1 to A2 is below 180 degrees so the Equation 1 of Dot product have produced it wrong
        theta_New = 2*(np.pi) -theta # Theta Clockwise. Because the dot product Eq.1 have produced the angle that is θ<180, and we desire the angle that is θ>=180


    #Up here we have found the angle θ of rotation to transform LCS into GCS.
    #Hypothetical scenario where the LCS and GCS are parallel.
    
    rotation_matrix = rotation_about_z_axis(theta_New)
    i=0
    Vert_Correct_GCS =  np.empty((0, Vert_Hyp_GCS.shape[1]))
    for i in range(8):
        Transformation2 = np.dot(rotation_matrix,Vert_Hyp_GCS[i])
        Vert_Correct_GCS = np.vstack((Vert_Correct_GCS, Transformation2))
        i = i+1
    print("The new matrix demonstrates the rotated IfcElement vertices from the Hypothetical scenario in GCS to the correct orientation regarding LCS",Vert_Correct_GCS )
    print("Step 4 Started: Translate the vertices back to LCS Origin point Oc  ")
    #Step 4: TRANSLATE BACK: 
    #   a) So far we have managed to rotate the IfcElement's vertices around z-axis in rotation point O(0,0,0) of GCS about rotation angle θ that specifies the 
    #       orientation of Local Coordinate system (LCS) of IfcElements
    #   b) Finally is time to translate the IfcElement back to the point Oc(3061,-1786,0) which is the origin point of IfcElement's LCS.
    #       Thus we transform the Matrix: Vert_Correct_GCS into Vert_Correct_LCS by adding the coordinates of the Oc(x,y,z)= Oc(3061,1786,0)
    #       Calculations: translation_matrix([[1., 0., 0., Dx],
                                                #[0., 1., 0., Dy],
                                                #[0., 0., 1., Dz],
                                                #[0., 0., 0., 1.]])*(Ax,Ay,Az,1)= (Ax + Dx, Ay + Dy, Az + Dz, 1)
    
    # Define translation vector
    translation = np.array([x, y, z]) #this command specifies the translation vector M2=[tx,ty,tz,1]
    # Create transformation matrix for translation
    translation_matrix[:3, 3] = translation #this command adds the 3x1 translation vector M2 in the 
    #translation matrix M, therefore gives values to the tx, ty,tz.
    #array([[1., 0., 0., Dx.],
    #[0., 1., 0., Dy.],
    #[0., 0., 1., Dz.],
    #[0., 0., 0., 1.]])
    # The Transformatino equation has certain structure as we explained before, thus, we need to transform point by point. This illustrated below.
    i=0
    Vert_Correct_LCS = np.empty((0, Vert_Correct_GCS.shape[1]))
    for i in range(8):
        Transformation3 = np.dot(translation_matrix,Vert_Correct_GCS[i] )
        Vert_Correct_LCS = np.vstack((Vert_Correct_LCS, Transformation3))
        i = i+1
    print("The new matrix demonstrates the translated IfcElement vertices from the Correct oriented in terms of GCS to correct oriented to LCS:",Vert_Correct_LCS )
    return Vert_Correct_LCS

patternIfcWall = r"IfcWall"
patternIfcColumn = r"IfcColumn"
patternIfcBeam = r"IfcBeam"
patternIfcSlab = r"IfcSlab"
patternIfcPipeSegment = r"IfcPipeSegment"
patternMappedItem = r"IfcMappedItem"
patternIfcExtrudedAreaSolid = r"IfcExtrudedAreaSolid"
patternIfcPolyline = r"IfcPolyline"
patternArbitary = r"IfcArbitraryClosedProfileDef"
patternRectangle = r"IfcRectangleProfileDef"
patternIfcOpeningElement = r"IfcOpeningElement"
patternClippingResultItem = r"IfcBooleanClippingResult"
slab = None
#rf: RefDirection, a: IfcArbitraryClosedProfileDef, r: IfcRectangleProfileDef, #w: IfcWalls, c: IfcColumns,etc
rfw = 0 
aw = 0
rw = 0
rfc = 0
ac = 0
rc = 0
rb = 0
N_Wall = 0
N_Column = 0
N_Beam = 0
N_Slab = 0
N_Pipe = 0
N_OpeningElement = 0
IfcWall_Matrices = []
IfcColumn_Matrices = []
IfcBeam_Matrices = []
IfcSlab_Matrices = []
IfcPipeSegment_Matrices = []
IfcOpeningElement_Matrices = []

faces = np.array([
                    [0, 1, 5, 4],  # 1st XZ face
                    [1, 2, 6, 5],  # 1st YZ face
                    [3, 0, 1, 2],  # 1st XY face
                    [2, 3, 7, 6],  # 2st XZ face
                    [3, 0, 4, 7],  # 2st YZ face
                    [7, 4, 5, 6]   # 2st XY face
                    ])
 # Create the Complete matrix that will contain the faces of all the elements
Complete_Matrix = np.array([["IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face 1 XZ", "Face 2 YZ", "Face 3 XY", "Face 4 XZ", "Face 5 YZ", "Face 6 XY"]], dtype=object)
print(Complete_Matrix)
Slab_Z_Location = np.array([["IfcType", "GUID", "Floor Level (z)"]], dtype=object)
for Object in IfcObjects:
            GUid = Object.GlobalId
            print(GUid)
            IfcIdentifier = Object.Tag
            #if GUid== "3XTMDtxRr2KPlqhaDaoB1B": #this is the rotated wall Floor 1
            #if object.GlobalId == "3XvVbJdofBkOIXzGLeWjIT": #this is parallel wall, but Horizontal
            if re.search(patternIfcWall, str(Object)):
                N_Wall = N_Wall + 1
                #if object.GlobalId == "257jRAh$nCtvOAJaAFLkbA": #this is the rotated wall 
                #if object.GlobalId == "3NU4ZQWO5798CgUAkHe54T": #this is the rotated wall 
                #if  re.search(patternIfcPipeSegment, str(object)):
                #if object.GlobalId == "24FIDMrNf7AR8O6xIs_9h7": #IfcColumn
                #if object.GlobalId == "3XvVbJdofBkOIXzGLeWjHX": #This a water network pipe with arbitary orientation
                #if slab.GlobalId == "2lWsAMV3j1TQjVIbH7IoN7": # This slab has Arbitary orientation and located in z=6096
                #if slab.GlobalId == "15dbPLPrf4WADu8KGE529J": # This slab has horizontal orientation and located in z=0
                #if slab.GlobalId == "2q1FgXqc9E6AIomE8nzLcR": # This slab has arbitary orientation and located in z=2900
                Local_Origin_IfcWall, Wall_RefDirection1, Wall_DirectionRatios, Wall_Axis1, Wall_Location2 = object_ObjectPlacement_inheritence(Object)
                #Local_Origin_IfcWall: This is the origin point of IfcWall's Local placement in terms of GCS 
                #wall_RefDirection1 this is the entity IfcDirection
                #Wall_DirectionRatios this is the attribute of entity IfcDirection
                Representations = object_Representations(Object)
                IfcGeometricRepresentationContext(Representations)
                print("The Object's Representations are:", Representations)
                
                Item_IfcExtrudedAreaSolid, Item_IfcIfcPolyline, PolylinePoints = ShapeRepresentations_loop(Representations)
                if Item_IfcExtrudedAreaSolid is not None:
                    SweptArea, Position, ExtrudedDirection, Depth, Location_of_IfcExtrudedAreaSolid, Axis, RefDirection = Item_IfcExtrudedAreaSolid_analysis(Item_IfcExtrudedAreaSolid)
                    # Here we distinguish between IfcElements with LCS Rotated in terms of GCS : IfcArbitraryClosedProfileDef, or parallel with GCS: IfcRectangleProfileDef
                    
                    #Step 1: a) Define Local origin point of Ifc Element in terms of GCS: Local_Orgin_Coord
                    #        b) Define the Vertices of the IfcElement in hypotherical scenario where the LCS is parallel with GCS 
                    print(Local_Origin_IfcWall)
                    #Here we store the Coordinates of the Local Origin point of the IfcWall in terms of the GCS
                    Local_Orgin_Coord = np.array(Local_Origin_IfcWall)
                    #To account for objects in different floors we adopt the cooridnates of IfcBuildingStorey's Origin point
                    #This point always refer to the (0,0,0) and in our cases differentiates only in z-coord, 
                    #So we inherit the z-coordinate of the floor to the origin point of its IfcObject
                    Abs_Coordinates = Wall_Location2.Coordinates # These are IFcBuildingStorey;s Origin point Coordinates.
                    print(f"Obs_Coordinates:",Abs_Coordinates)
                    Abs_Coordinates = np.array(Abs_Coordinates)
                    Floor_z_Coord = Abs_Coordinates[2] #We take the Floor level of the current IfcObject by its parent object (IfcBuildingStorey)
                    print(Floor_z_Coord)
                    x = Local_Orgin_Coord[0]
                    y = Local_Orgin_Coord[1]
                    z = Local_Orgin_Coord[2] # instead of Floor_z_Coord. Floor_z_Coord creates an offset of the value of Local_Orgin_Coord[2], becasue if the z coordinate of the wall is non zero then the vlues of z coordinate represents the height offset of the wall in terms of the Building storey level.
                    Floor_Z_Coord = Local_Orgin_Coord[2]
                    if Wall_RefDirection1 is None:
                        rfw = rfw +1
                        print("The IfcElement has LCS parallel with GCS so there is no need for transformation of its vertices ")
                        Wall_XDim, Wall_YDim = Item_SweptArea_Rectangle(SweptArea) #XDim, and YDim, define the profile area of the wall
                        print(Wall_XDim, Wall_YDim)
                        #Find the Coordinates of the IfcElement vertices in terms of GCS in the Hypothetical scenario where LGC is parallel with GCS
                        Ax = x
                        Ay = y - Wall_YDim/2
                        Bx = x + Wall_XDim
                        By = Ay
                        Cx = Bx 
                        Cy = y + Wall_YDim/2
                        Dx = x 
                        Dy = Cy
                        Vert_Correct_LCS = [(Ax, Ay, z,1),
                                        (Bx, By, z,1),
                                        (Cx, Cy, z,1),
                                        (Dx, Dy, z,1),
                                        (Ax, Ay, z + Depth,1),
                                        (Bx, By, z + Depth,1),
                                        (Cx, Cy, z + Depth,1),
                                        (Dx, Dy, z + Depth,1)]
                        Vert_Correct_LCS = np.array(Vert_Correct_LCS)
                        #Then if RefDirection is None the ifc element might be IfcSlab, then in this case:
                    else:
                        if re.search(patternArbitary, str(SweptArea)):
                            aw = aw + 1
                            print("Step 1: Find the IfcElement Vertices in the hypothetical scenario in point Oc where the LGC is parallel with GCS.")
                            #For IfcWall we have 8 vertices. Here we define the vertices as if the column's LCS were in parallel with the GCS.
                            #the nomenclature goes like this: A symbolises the name of the vertex and 1 symbolises the level in which the vertex exists, 
                            #we will always have 2 levels ()
                            LGC_direction_Vect = Wall_DirectionRatios
                            ProfDefPoints = Item_SweptArea_Arbitary(SweptArea)
                            Point_A = ProfDefPoints[0].Coordinates
                            Point_A = np.array(Point_A)
                            print("Coordinates of point A of IFCWall area:", Point_A)
                            Point_B = ProfDefPoints[1].Coordinates
                            Point_B = np.array(Point_B)
                            print("Coordinates of point B of IFCWall area:", Point_B)
                            Point_C = ProfDefPoints[2].Coordinates
                            Point_C = np.array(Point_C)
                            print("Coordinates of point C of IFCWall area:", Point_C)
                            Point_D = ProfDefPoints[3].Coordinates
                            Point_D = np.array(Point_D)
                            print("Coordinates of point D of IFCWall area:", Point_D)
                            #find absolute difference of Ay and By coordinate, it will be used later in the code
                            Abs_Dif_ABy = abs(Point_B[1]-Point_A[1]) 
                            print("Absolute difference, IfcWall thickness:", Abs_Dif_ABy)
                            #the IfcElement's vertices in terms of GCS in hypothetical scenario are (26/03/2024)
                            Vert_Hyp_LCS = [(x+ Point_A[0], y + Point_A[1], z,1),
                                            (x+ Point_B[0], y + Point_B[1], z,1),
                                            (x+ Point_C[0], y + Point_C[1], z,1),
                                            (x+ Point_D[0], y + Point_D[1], z,1),
                                            (x+ Point_A[0], y + Point_A[1], z + Depth,1),
                                            (x+ Point_B[0], y + Point_B[1], z + Depth,1),
                                            (x+ Point_C[0], y + Point_C[1], z + Depth,1),
                                            (x+ Point_D[0], y + Point_D[1], z + Depth,1)] #The coordinate z is derived from IfcLocalPlacement of the IfcElement
                            Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                            print(Vert_Hyp_LCS)
                            Vert_Correct_LCS = Step2_3_4_Vertices_Transformation(Vert_Hyp_LCS, x, y , z, LGC_direction_Vect)
                        if re.search(patternRectangle, str(SweptArea)):
                            rw = rw + 1
                            print("Step 1: Find the IfcElement Vertices in the hypothetical scenario in point Oc where the LGC is parallel with GGC.")
                            #Step1: For IfcColumn we have 8 vertices. Here we define the vertices as if the column's LCS were in parallel with the GCS.
                            #the nomenclature goes like this: A symbolises the name of the vertex and 1 symbolises the level in which the vertex exists, 
                            #we will always have 2 levels ()
                            LGC_direction_Vect = Wall_DirectionRatios #This defines the vertex A2 (RefDirection), namely the direction of the Local positive X-axis of Ifcwall
                            Wall_XDim, Wall_YDim = Item_SweptArea_Rectangle(SweptArea) #XDim, and YDim, define the profile area of the wall
                            print(Wall_XDim, Wall_YDim)
                            #Find the Coordinations of the IfcElement vertices in terms of GCS in the Hypothetical scenario where LGC is parallel with GCS
                            #And vertices are located in origin point Oc.
                            Ax = x
                            Ay = y - Wall_YDim/2
                            Bx = x + Wall_XDim
                            By = Ay
                            Cx = Bx 
                            Cy = y + Wall_YDim/2
                            Dx = x 
                            Dy = Cy
                            Vert_Hyp_LCS = [(Ax, Ay, z,1),
                                            (Bx, By, z,1),
                                            (Cx, Cy, z,1),
                                            (Dx, Dy, z,1),
                                            (Ax, Ay, z + Depth,1),
                                            (Bx, By, z + Depth,1),
                                            (Cx, Cy, z + Depth,1),
                                            (Dx, Dy, z + Depth,1)]
                            Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                            print("The Vertces in Step 1 are:", Vert_Hyp_LCS)
                            Vert_Correct_LCS = Step2_3_4_Vertices_Transformation(Vert_Hyp_LCS, x, y , z, LGC_direction_Vect)
                    
                    # Delete the last column from all rows
                    Vert_Correct_LCS = np.delete(Vert_Correct_LCS, -1, axis=1)
                    print(Vert_Correct_LCS)
                    # Creating a list of the matrices, its GUID, and its ID
                    IfcWall_Matrices.append((Vert_Correct_LCS, GUid, N_Wall, z))
                else:
                    print("this wall is not represented by SweptSolid Representation Type, so it is excluded from the analysis")
                if N_Wall == 1: #Create those two variables only one time, and use them for every IfcObject.
                    faces = np.array([
                    [0, 1, 5, 4],  # 1st XZ face
                    [1, 2, 6, 5],  # 1st YZ face
                    [3, 0, 1, 2],  # 1st XY face
                    [2, 3, 7, 6],  # 2st XZ face
                    [3, 0, 4, 7],  # 2st YZ face
                    [7, 4, 5, 6]   # 2st XY face
                    ])
                    triangles = np.array([
                        [0, 1, 5 ],  # Tr1 1st XZ face
                        [5, 4, 0 ],  # Tr2 1st XZ face 
                        [1, 2, 6 ],  # Tr1 1st YZ face
                        [6, 5, 1 ],  # Tr2 1st YZ face
                        [0, 1, 2 ],  # Tr1 1st XY face
                        [2, 3, 0 ],  # Tr2 1st XY face
                        [2, 3, 7 ],  # Tr1 2nd XZ face #Second face starts
                        [7, 6, 2 ],  # Tr2 2nd XZ face 
                        [3, 0, 4 ],  # Tr1 2nd YZ face
                        [4, 7, 3 ],  # Tr2 2nd YZ face
                        [7, 4, 5 ],  # Tr1 2nd XY face
                        [5, 6, 7 ],  # Tr2 2nd XY face
                    ])
                # Separate faces into separate arrays
                separate_faces = [Vert_Correct_LCS[face_indices] for face_indices in faces]
                New_object = np.array([patternIfcWall, GUid, IfcIdentifier, N_Wall, z , None, None, None, None, None, None], dtype=object)
                #Iterate through the number of faces 
                print(separate_faces)
                for i, face in enumerate(separate_faces):    
                    New_object[i + 5] = face  # Assign the first value of each array to the corresponding column in new_row
                # Append the new row to the array
                Complete_Matrix = np.vstack([Complete_Matrix, New_object])
                print(Complete_Matrix)

                #For IfcWall's Opening Elements
                voids = ifc_file.get_inverse(Object, "HasOpenings")
                print(voids)
                RelVoidsElements = []
                for void in voids:
                    if void.is_a("IfcRelVoidsElement"):
                        RelVoidsElements.append(void)
                print(RelVoidsElements) #now we have imported in this array all the IfcRelVoidsElement that are related to the current IfcWall
                #A IfcRelVoidsElement is associated with the openning element through the relationship "RelatedOpeningElement" which leads to IfcOpeningElement class
                OpeningElements = []
                for RelVoidsElement in RelVoidsElements:
                    OpeningElements.append(RelVoidsElement.RelatedOpeningElement)
                print(OpeningElements) # Now we have the three desirable IfcOpeningElements that are associated with the current wall, so we can include them in the transformaton process later on
                #Now we can start parsing each IfcOpeningElement, we only need its IfcLocalplacement to find the relation with the wall. As we argued the IfcOpenning Element
                #is related to IfcWall's LocalPlacement regarding IFcDocumentation. So we only need the first IfcLocalplacement of the IfcOpeningElement.
                #Then we continue the process by finding the vertices of IFCWall and the vertices of each of IfcOpeningElement, we use the transformation function seperately for each entity...
                #Check the two different cases where refdirection of the wall is positive and negative
                N_OpeningElement = N_Wall
                for OpeningElement in OpeningElements:
                        N_OpeningElement = round(N_OpeningElement + 0.1, 1) #Round the ID in one decimal point
                        GUid_OE = OpeningElement.GlobalId
                        print(GUid)
                    #if GUid_OE == "05hIZEszH5GBM5fTU338u2":
                    #if GUid_OE == "05hIZEszH5GBM5fSI338v_":
                        IfcIdentifier_OE = OpeningElement.Tag
                        Local_Origin_IfcOpeningElement, OE_RefDirection1, OE_DirectionRatios, Axis1, Location2 = object_ObjectPlacement_inheritence(OpeningElement)
                        #from above results we need only the Local_Origin_Ifcobject variable of OpeningElement
                        # Regarding IfcOpeningElement entity, there is always one IfcShapeRepresentation and is always IfcExtrudedAreaSolid with IfcRectangleProfileDef
                        Representations = object_Representations(OpeningElement)
                        IfcGeometricRepresentationContext(Representations)
                        print("The Object's Representations:", Representations)
                        
                        Item_IfcExtrudedAreaSolid, Item_IfcIfcPolyline, PolylinePoints = ShapeRepresentations_loop(Representations)
                        if Item_IfcExtrudedAreaSolid is not None:
                            SweptArea, Position, ExtrudedDirection, Depth, Location_of_IfcExtrudedAreaSolid, Axis, RefDirection = Item_IfcExtrudedAreaSolid_analysis(Item_IfcExtrudedAreaSolid)
                            # Here we distinguish between IfcElements with LCS Rotated in terms of GCS : IfcArbitraryClosedProfileDef, or parallel with GCS: IfcRectangleProfileDef
                            
                            #Step 1: a) Define Local origin point of Ifc Element in terms of GCS: Local_Orgin_Coord
                            #        b) Define the Vertices of the IfcElement in hypotherical scenario where the LCS is parallel with GCS 
                            print(Local_Origin_IfcOpeningElement)
                            #Here we store the Coordinates of the Local Origin point of the IfcWall in terms of the GCS
                            Local_Orgin_Coord_IfcWall= np.array(Local_Origin_IfcWall)
                            x_wall = Local_Orgin_Coord_IfcWall[0]
                            y_wall = Local_Orgin_Coord_IfcWall[1]
                            z_wall = Local_Orgin_Coord_IfcWall[2]
                            #Here we store the Coordinates of the Local Origin point of the IfcOpeningElement in terms of the ifcWall's LCS
                            Local_Orgin_Coord_IfcOpeningElement= np.array(Local_Origin_IfcOpeningElement)
                            x_OE = Local_Orgin_Coord_IfcOpeningElement[0]
                            y_OE= Local_Orgin_Coord_IfcOpeningElement[1]
                            XDim, YDim = Item_SweptArea_Rectangle(SweptArea) #XDim, and YDim, define the profile area of the wall
                            print(XDim, YDim)
                            #We distinguish two Cases: Case 1: The IfcWall's X-axis (RefDirection attribute) is located in one of the GCS positive axis X or Y
                            #Case 2: The IfcWall's X-axis (RefDirection attribute) is located in one of the GCS negative axis X or Y
                            IfcWall_RefDirection_Coordinates = np.array(Wall_DirectionRatios) 
                            IfcOE_RefDirection_Coordinates = np.array(OE_DirectionRatios) 
                            if Wall_RefDirection1 is None or np.all(IfcWall_RefDirection_Coordinates >=  0):
                                if  OE_RefDirection1 is None or np.all(IfcOE_RefDirection_Coordinates >=  0):
                                    #The Opening elemnt'svertices are finding just like in case 1.2
                                    Ax = x_wall + x_OE
                                    Ay = y_wall + y_OE
                                    Bx = Ax + YDim
                                    By = Ay
                                    Cx = Bx 
                                    Cy = By + Wall_YDim
                                    Dx = Ax
                                    Dy = Cy
                                    Vert_Hyp_LCS = [(Ax, Ay, z_wall,1),
                                                    (Bx, By, z_wall,1),
                                                    (Cx, Cy, z_wall,1),
                                                    (Dx, Dy, z_wall,1),
                                                    (Ax, Ay, z_wall + XDim,1),
                                                    (Bx, By, z_wall + XDim,1),
                                                    (Cx, Cy, z_wall + XDim,1),
                                                    (Dx, Dy, z_wall + XDim,1)]
                                    Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                                elif np.any(IfcOE_RefDirection_Coordinates < 0):
                                    #Case 1.1
                                    Dx = x_wall + x_OE
                                    Dy = y_wall + y_OE
                                    Cx = Dx + YDim
                                    Cy = Dy
                                    Bx = Cx 
                                    By = Dy - Wall_YDim
                                    Ax = Dx
                                    Ay = By
                                    Vert_Hyp_LCS = [(Ax, Ay, z_wall,1),
                                                    (Bx, By, z_wall,1),
                                                    (Cx, Cy, z_wall,1),
                                                    (Dx, Dy, z_wall,1),
                                                    (Ax, Ay, z_wall + XDim,1),
                                                    (Bx, By, z_wall + XDim,1),
                                                    (Cx, Cy, z_wall + XDim,1),
                                                    (Dx, Dy, z_wall + XDim,1)]
                                    Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                            else: #Case 2.1
                                if OE_RefDirection1 is None: #the case is similar to the Case 1.2
                                    Ax = x_wall + x_OE
                                    Ay = y_wall + y_OE
                                    Bx = Ax + YDim
                                    By = Ay
                                    Cx = Bx 
                                    Cy = By + Wall_YDim
                                    Dx = Ax
                                    Dy = Cy
                                    Vert_Hyp_LCS = [(Ax, Ay, z_wall,1),
                                                    (Bx, By, z_wall,1),
                                                    (Cx, Cy, z_wall,1),
                                                    (Dx, Dy, z_wall,1),
                                                    (Ax, Ay, z_wall + XDim,1),
                                                    (Bx, By, z_wall + XDim,1),
                                                    (Cx, Cy, z_wall + XDim,1),
                                                    (Dx, Dy, z_wall + XDim,1)]
                                    Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                                elif np.all(IfcOE_RefDirection_Coordinates >= 0): # Case 2.2
                                #Find the Coordinates of the IfcElement vertices in terms of GCS in the Hypothetical scenario where LGC is parallel with GCS
                                #And vertices are located in origin point Oc.
                                    Bx = x_wall + x_OE
                                    By = y_wall + y_OE
                                    Cx = Bx
                                    Cy = By + Wall_YDim
                                    Dx = Cx - YDim
                                    Dy = Cy
                                    Ax = Dx
                                    Ay = By
                                    Vert_Hyp_LCS = [(Ax, Ay, z_wall,1),
                                                    (Bx, By, z_wall,1),
                                                    (Cx, Cy, z_wall,1),
                                                    (Dx, Dy, z_wall,1),
                                                    (Ax, Ay, z_wall + XDim,1),
                                                    (Bx, By, z_wall + XDim,1),
                                                    (Cx, Cy, z_wall + XDim,1),
                                                    (Dx, Dy, z_wall + XDim,1)]
                                    Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                                    
                                else: #Case 2.3
                                    Cx = x_wall + x_OE
                                    Cy = y_wall + y_OE
                                    Dx = Cx - YDim
                                    Dy = Cy
                                    Ax = Dx 
                                    Ay = Dy - Wall_YDim
                                    Bx = Cx
                                    By = Ay
                                    Vert_Hyp_LCS = [(Ax, Ay, z_wall,1),
                                                    (Bx, By, z_wall,1),
                                                    (Cx, Cy, z_wall,1),
                                                    (Dx, Dy, z_wall,1),
                                                    (Ax, Ay, z_wall + XDim,1),
                                                    (Bx, By, z_wall + XDim,1),
                                                    (Cx, Cy, z_wall + XDim,1),
                                                    (Dx, Dy, z_wall + XDim,1)]
                                    Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                            if Wall_RefDirection1 is None:
                                print("The IfcElement has LCS parallel with GCS so there is no need for transformation of its vertices ")
                                #Then if RefDirection is None the ifc element might be IfcSlab, then in this case:
                                Vert_Correct_LCS = Vert_Hyp_LCS
                            else: 
                                print("Step 1: Find the IfcElement Vertices in the hypothetical scenario in point Oc where the LGC is parallel with GGC.")
                                #Step1: For IfcColumn we have 8 vertices. Here we define the vertices as if the column's LCS were in parallel with the GCS.
                                #the nomenclature goes like this: A symbolises the name of the vertex and 1 symbolises the level in which the vertex exists, 
                                #we will always have 2 levels ()
                                LGC_direction_Vect = Wall_DirectionRatios #This defines the vertex A2 (RefDirection), namely the direction of the Local positive X-axis of Ifcwall
                                print("The Vertces in Step 1 are:", Vert_Hyp_LCS)
                                Vert_Correct_LCS = Step2_3_4_Vertices_Transformation(Vert_Hyp_LCS, x_wall, y_wall , z_wall, LGC_direction_Vect)
                            
                            # Delete the last column from all rows
                            Vert_Correct_LCS = np.delete(Vert_Correct_LCS, -1, axis=1)
                            print(Vert_Correct_LCS)
                            # Creating a list of the matrices, its GUID, and its ID
                            IfcOpeningElement_Matrices.append((Vert_Correct_LCS, GUid_OE, IfcIdentifier_OE, N_OpeningElement, z_wall))
                        else:
                            print("this OpeningElement is not represented by SweptSolid Representation Type, so it is excluded from the analysis")
                        # Separate faces into separate arrays
                        separate_faces = [Vert_Correct_LCS[face_indices] for face_indices in faces]
                        New_object = np.array([patternIfcOpeningElement, GUid_OE, IfcIdentifier_OE, N_OpeningElement, z_wall , None, None, None, None, None, None], dtype=object)
                        #Iterate through the number of faces 
                        print(separate_faces)
                        for i, face in enumerate(separate_faces):    
                            New_object[i + 5] = face  # Assign the first value of each array to the corresponding column in new_row
                        # Append the new row to the array
                        Complete_Matrix = np.vstack([Complete_Matrix, New_object])
                        print(Complete_Matrix)
            elif  re.search(patternIfcColumn, str(Object)):
                N_Column = N_Column + 1
                Local_Origin_Ifcobject_coords, RefDirection1, DirectionRatios, Axis1, Location2 = object_ObjectPlacement_inheritence(Object)
                Representations = object_Representations(Object)
                IfcGeometricRepresentationContext(Representations)
                print("The Object's Representations:", Representations)
                
                Item_IfcExtrudedAreaSolid, Item_IfcIfcPolyline, PolylinePoints = ShapeRepresentations_loop(Representations)
                if Item_IfcExtrudedAreaSolid is not None:
                    SweptArea, Position, ExtrudedDirection, Depth, Location_of_IfcExtrudedAreaSolid, Axis_of_IfcExtrudedAreaSolid, RefDirection_of_IfcExtrudedAreaSolid = Item_IfcExtrudedAreaSolid_analysis(Item_IfcExtrudedAreaSolid)
                    # Here we distinguish between IfcElements with LCS Rotated in terms of GCS : IfcArbitraryClosedProfileDef, or parallel with GCS: IfcRectangleProfileDef
                    
                    #Step 1: a) Define Local origin point of Ifc Element in terms of GCS: Local_Orgin_Coord
                    #        b) Define the Vertices of the IfcElement in hypotherical scenario where the LCS is parallel with GCS 
                    print(Local_Origin_Ifcobject_coords)
                    print(Axis1)
                    print(RefDirection1)
                    Local_Origin_Ifcobject_coords = np.array(Local_Origin_Ifcobject_coords)
                    if Axis1 == None and RefDirection1 == None and np.all(Local_Origin_Ifcobject_coords) == 0: # This is a rarely case in which the Object placement of the Column is specified by the local placement of IFcExtrudedAreasSolid
                        Location_of__coords_IfcExtrudedAreaSolid = Location_of_IfcExtrudedAreaSolid.Coordinates
                        Location_of__coords_IfcExtrudedAreaSolid = np.array(Location_of__coords_IfcExtrudedAreaSolid)
                        x = Location_of__coords_IfcExtrudedAreaSolid[0]
                        y = Location_of__coords_IfcExtrudedAreaSolid[1]
                        z = Location_of__coords_IfcExtrudedAreaSolid[2]
                    else:
                        #Here we store the Coordinates of the Local Origin point of the IfcWall in terms of the GCS
                        
                        #To account for objects in different floors we adopt the cooridnates of IfcBuildingStorey's Origin point
                        #This point always refer to the (0,0,0) and in our cases differentiates only in z-coord, 
                        #So we inherit the z-coordinate of the floor to the origin point of its IfcObject
                        Obs_Coordinates = Location2.Coordinates # These are IFcBuildingStorey;s Origin point Coordinates.
                        print(f"Obs_Coordinates:",Obs_Coordinates)
                        Obs_Coordinates = np.array(Obs_Coordinates)
                        Floor_z_Coord = Obs_Coordinates[2]
                        print(Floor_z_Coord)
                        x = Local_Origin_Ifcobject_coords[0]
                        y = Local_Origin_Ifcobject_coords[1]
                        z = Floor_z_Coord
                    
                    if RefDirection1 is None:
                        rfc= rfc+1
                        print("The IfcElement has LCS parallel with GCS so there is no need for transformation of its vertices ")
                        XDim, YDim = Item_SweptArea_Rectangle(SweptArea) #XDim, and YDim, define the profile area of the wall
                        print(XDim, YDim)
                        #Find the Coordinations of the IfcElement vertices in terms of GCS in the Hypothetical scenario where LGC is parallel with GCS
                        Ax = x + YDim/2
                        Ay = y - XDim/2
                        Bx = x + YDim/2
                        By = y + XDim/2
                        Cx = x - YDim/2
                        Cy = y + XDim/2
                        Dx = Cx
                        Dy = Ay
                        Vert_Correct_LCS = [(Ax, Ay, z,1), #Here the correct vertices coordinates are found directly.
                                            (Bx, By, z,1),
                                            (Cx, Cy, z,1),
                                            (Dx, Dy, z,1),
                                            (Ax, Ay, z + Depth,1),
                                            (Bx, By, z + Depth,1),
                                            (Cx, Cy, z + Depth,1),
                                            (Dx, Dy, z + Depth,1)]
                        Vert_Correct_LCS = np.array(Vert_Correct_LCS)
                    else:
                        if re.search(patternArbitary, str(SweptArea)):
                            ac = ac+1
                            print("Step 1: Find the IfcElement Vertices in the hypothetical scenario in point Oc where the LGC is parallel with GCS.")
                            #For IfcWall we have 8 vertices. Here we define the vertices as if the column's LCS were in parallel with the GCS.
                            #the nomenclature goes like this: A symbolises the name of the vertex and 1 symbolises the level in which the vertex exists, 
                            #we will always have 2 levels ()
                            LGC_direction_Vect = DirectionRatios
                            ProfDefPoints = Item_SweptArea_Arbitary(SweptArea)
                            Point_A = ProfDefPoints[0].Coordinates
                            Point_A = np.array(Point_A)
                            print("Coordinates of first point A of IFCWall area:", Point_A)
                            Point_B = ProfDefPoints[1].Coordinates
                            Point_B = np.array(Point_B)
                            print("Coordinates of first point B of IFCWall area:", Point_B)
                            Point_C = ProfDefPoints[2].Coordinates
                            Point_C = np.array(Point_C)
                            print("Coordinates of first point B of IFCWall area:", Point_C)
                            Point_D = ProfDefPoints[3].Coordinates
                            Point_D = np.array(Point_D)
                            print("Coordinates of first point B of IFCWall area:", Point_D)
                            #find absolute difference of Ay and By coordinate, it will be used later in the code
                            Abs_Dif_ABy = abs(Point_B[1]-Point_A[1]) 
                            print("Absolute difference, IfcWall thickness:", Abs_Dif_ABy)
                            #the IfcElement's vertices in terms of GCS in hypothetical scenario are (26/03/2024)
                            Vert_Hyp_LCS = [(x+ Point_A[0], y + Point_A[1], z,1),
                                            (x+ Point_B[0], y + Point_B[1], z,1),
                                            (x+ Point_C[0], y + Point_C[1], z,1),
                                            (x+ Point_D[0], y + Point_D[1], z,1),
                                            (x+ Point_A[0], y + Point_A[1], z + Depth,1),
                                            (x+ Point_B[0], y + Point_B[1], z + Depth,1),
                                            (x+ Point_C[0], y + Point_C[1], z + Depth,1),
                                            (x+ Point_D[0], y + Point_D[1], z + Depth,1)] #The coordinate z is derived from IfcLocalPlacement of the IfcElement
                            Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                            print(Vert_Hyp_LCS)
                            Vert_Correct_LCS = Step2_3_4_Vertices_Transformation(Vert_Hyp_LCS, x, y , z, LGC_direction_Vect)
                        if re.search(patternRectangle, str(SweptArea)):
                            rc = rc+1
                            print("Step 1: Find the IfcElement Vertices in the hypothetical scenario in point Oc where the LGC is parallel with GGC.")
                            #Step1: For IfcColumn we have 8 vertices. Here we define the vertices as if the column's LCS were in parallel with the GCS.
                            #the nomenclature goes like this: A symbolises the name of the vertex and 1 symbolises the level in which the vertex exists, 
                            #we will always have 2 levels ()
                            LGC_direction_Vect = DirectionRatios #This defines the vertex A2 (RefDirection), namely the direction of the Local positive X-axis of Ifcwall
                            XDim, YDim = Item_SweptArea_Rectangle(SweptArea) #XDim, and YDim, define the profile area of the wall
                            print(XDim, YDim)
                            #Find the Coordinations of the IfcElement vertices in terms of GCS in the Hypothetical scenario where LGC is parallel with GCS
                            #And vertices are located in origin point Oc.
                            Ax = x + YDim/2
                            Ay = y - XDim/2
                            Bx = x + YDim/2
                            By = y + XDim/2
                            Cx = x - YDim/2
                            Cy = y + XDim/2
                            Dx = Cx
                            Dy = Ay
                            Vert_Hyp_LCS = [(Ax, Ay, z,1),
                                            (Bx, By, z,1),
                                            (Cx, Cy, z,1),
                                            (Dx, Dy, z,1),
                                            (Ax, Ay, z + Depth,1),
                                            (Bx, By, z + Depth,1),
                                            (Cx, Cy, z + Depth,1),
                                            (Dx, Dy, z + Depth,1)]
                            Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                            print("The Vertces in Step 1 are:", Vert_Hyp_LCS)
                            Vert_Correct_LCS = Step2_3_4_Vertices_Transformation(Vert_Hyp_LCS, x, y ,z , LGC_direction_Vect)
                    # Delete the last column from all rows
                    Vert_Correct_LCS = np.delete(Vert_Correct_LCS, -1, axis=1)
                    print(Vert_Correct_LCS)
                    # Creating a list of the matrices, its GUID, and its ID
                    IfcColumn_Matrices.append((Vert_Correct_LCS, GUid, N_Column, z ))
                else:
                    print("this Column is not represented by SweptSolid Representation Type, so it is excluded from the analysis")
                if N_Column == 1:
                    faces = np.array([
                    [0, 1, 5, 4],  # 1st XZ face
                    [1, 2, 6, 5],  # 1st YZ face
                    [3, 0, 1, 2],  # 1st XY face
                    [2, 3, 7, 6],  # 2st XZ face
                    [3, 0, 4, 7],  # 2st YZ face
                    [7, 4, 5, 6]   # 2st XY face
                    ])
                # Separate faces into separate arrays
                separate_faces = [Vert_Correct_LCS[face_indices] for face_indices in faces]
                New_object = np.array([patternIfcColumn, GUid, IfcIdentifier, N_Column, z, None, None, None, None, None, None], dtype=object)
                #Iterate through the number of faces 
                for i, face in enumerate(separate_faces):
                    if i==0: #face 1 XZ in Complete matrix
                        New_object[6] = face  # Assign the first value of each array to the corresponding column in new_rowNew_object[i + 5] = face  # Assign the first value of each array to the corresponding column in new_row
                    elif i==1: #Face 2 YZ in Complete matrix
                        New_object[5] = face  # Assign the first value of each array to the corresponding column in new_row
                    elif i==2: #face 3 XY in Complete matrix
                        New_object[7] = face
                    elif i==3: #face 4 XZ in Complete matrix
                        New_object[9] = face
                    elif i==4: #face 5 YZ in Complete matrix
                        New_object[8] = face
                    elif i==5: #face 6 XY in Complete matrix
                        New_object[10] = face
                
                # Append the new row to the array
                Complete_Matrix = np.vstack([Complete_Matrix, New_object])
                print(Complete_Matrix)           
            elif  re.search(patternIfcBeam, str(Object)): 
                N_Beam = N_Beam + 1
                Local_Origin_Ifcobject, RefDirection1, DirectionRatios, Axis1, Location2 = object_ObjectPlacement_inheritence(Object)
                Representations = object_Representations(Object)
                IfcGeometricRepresentationContext(Representations)
                print("The Object's Representations:", Representations)
                
                Item_IfcExtrudedAreaSolid, Item_IfcIfcPolyline, PolylinePoints = ShapeRepresentations_loop(Representations)
                if Item_IfcExtrudedAreaSolid is not None:
                    SweptArea, Position, ExtrudedDirection, Depth, Location_of_IfcExtrudedAreaSolid, Axis, RefDirection = Item_IfcExtrudedAreaSolid_analysis(Item_IfcExtrudedAreaSolid)
                    # Here we distinguish between IfcElements with LCS Rotated in terms of GCS : IfcArbitraryClosedProfileDef, or parallel with GCS: IfcRectangleProfileDef
                    
                    #Step 1: a) Define Local origin point of Ifc Element in terms of GCS: Local_Orgin_Coord
                    #        b) Define the Vertices of the IfcElement in hypotherical scenario where the LCS is parallel with GCS 
                    print(Local_Origin_Ifcobject)
                    #Here we store the Coordinates of the Local Origin point of the IfcWall in terms of the GCS
                    Local_Orgin_Coord = np.array(Local_Origin_Ifcobject)
                    
                    x = Local_Orgin_Coord[0]
                    y = Local_Orgin_Coord[1]
                    z = Local_Orgin_Coord[2]
                    #IfcBeam has only IfcRectangleprofileDef, so we exclude all other cases (RefDirection=None, and ArbitaryProf)
                    if re.search(patternRectangle, str(SweptArea)):           
                        rb = rb +1
                        print("Step 1: Find the IfcElement Vertices in the hypothetical scenario in point Oc where the LGC is parallel with GGC.")
                        #Step1: For IfcColumn we have 8 vertices. Here we define the vertices as if the column's LCS were in parallel with the GCS.
                        #the nomenclature goes like this: A symbolises the name of the vertex and 1 symbolises the level in which the vertex exists, 
                        #we will always have 2 levels ()
                        LGC_direction_Vect = DirectionRatios #This defines the vertex A2 (RefDirection), namely the direction of the Local positive X-axis of Ifcwall
                        XDim, YDim = Item_SweptArea_Rectangle(SweptArea) #XDim, and YDim, define the profile area of the wall
                        print(XDim, YDim)
                        Coordinates_of_FloorLevel = Location2.Coordinates
                        print(f"Coordinates_of_FloorLevel", Coordinates_of_FloorLevel)
                        Coordinates_of_FloorLevel = np.array(Coordinates_of_FloorLevel)
                        Floor_Z_Coord = Coordinates_of_FloorLevel[2]
                        print(f"Floor_Z_Coord", Floor_Z_Coord)
                        #Find the Coordinations of the IfcElement vertices in terms of GCS in the Hypothetical scenario where LGC is parallel with GCS
                        #And vertices are located in origin point Oc.
                        Ax = x + Depth/2
                        Ay = y - YDim/2
                        Bx = x + Depth/2
                        By = y + YDim/2
                        Cx = x - Depth/2
                        Cy = y + YDim/2
                        Dx = Cx
                        Dy = Ay
                        #We assume that the origin point of the beam goes through the mass center of the beam. and the mass center of its cross section
                        z1 = Floor_Z_Coord + z + XDim/2 # Here we find the vertices A1, B1, C1, D1 in the upper XY plane of the Beam
                        z2 = Floor_Z_Coord + z - XDim/2 # Here we find the vertices A2, B2, C2, D2 in the down XY plane of the Beam
                        Vert_Hyp_LCS = [(Ax, Ay, z1,1),
                                            (Bx, By, z1,1),
                                            (Cx, Cy, z1,1),
                                            (Dx, Dy, z1,1),
                                            (Ax, Ay, z2,1),
                                            (Bx, By, z2,1),
                                            (Cx, Cy, z2,1),
                                            (Dx, Dy, z2,1)]
                        Vert_Hyp_LCS = np.array(Vert_Hyp_LCS)
                        print("The Vertces in Step 1 are:", Vert_Hyp_LCS)
                        Vert_Correct_LCS = Step2_3_4_Vertices_Transformation(Vert_Hyp_LCS, x, y , z, LGC_direction_Vect)
                        # Delete the last column from all rows
                    Vert_Correct_LCS = np.delete(Vert_Correct_LCS, -1, axis=1)
                    print(Vert_Correct_LCS)
                    # Creating a list of the matrices, its GUID, and its ID
                    IfcBeam_Matrices.append((Vert_Correct_LCS, GUid, N_Beam, z1 ))
                else:
                    print("this Beam is not represented by SweptSolid Representation Type, so it is excluded from the analysis")                
                if N_Beam == 1:
                    faces = np.array([
                    [0, 1, 5, 4],  # 1st XZ face
                    [1, 2, 6, 5],  # 1st YZ face
                    [3, 0, 1, 2],  # 1st XY face
                    [2, 3, 7, 6],  # 2st XZ face
                    [3, 0, 4, 7],  # 2st YZ face
                    [7, 4, 5, 6]   # 2st XY face
                    ])
                # Separate faces into separate arrays
                separate_faces = [Vert_Correct_LCS[face_indices] for face_indices in faces]
                New_object = np.array([patternIfcBeam, GUid, IfcIdentifier, N_Beam, Floor_Z_Coord, None, None, None, None, None, None], dtype=object)
                #Iterate through the number of faces 
                for i, face in enumerate(separate_faces):
                    if i==0: #face 1 XZ in Complete matrix
                        New_object[6] = face  # Assign the first value of each array to the corresponding column in new_rowNew_object[i + 5] = face  # Assign the first value of each array to the corresponding column in new_row
                    elif i==1: #Face 2 YZ in Complete matrix
                        New_object[5] = face  # Assign the first value of each array to the corresponding column in new_row
                    elif i==2: #face 3 XY in Complete matrix
                        New_object[7] = face
                    elif i==3: #face 4 XZ in Complete matrix
                        New_object[9] = face
                    elif i==4: #face 5 YZ in Complete matrix
                        New_object[8] = face
                    elif i==5: #face 6 XY in Complete matrix
                        New_object[10] = face
                
                # Append the new row to the array
                Complete_Matrix = np.vstack([Complete_Matrix, New_object])
                print(Complete_Matrix)      
            elif  re.search(patternIfcSlab, str(Object)):  
                N_Slab = N_Slab + 1
                Local_Origin_Ifcobject, RefDirection1, DirectionRatios, Axis1, Location2 = object_ObjectPlacement_inheritence(Object)
                Representations = object_Representations(Object)
                IfcGeometricRepresentationContext(Representations)
                print("The Object's Representations:", Representations)
                
                Item_IfcExtrudedAreaSolid, Item_IfcIfcPolyline, SlabPoints = ShapeRepresentations_loop(Representations)
                if Item_IfcExtrudedAreaSolid is not None:
                    SweptArea, Position, ExtrudedDirection, Depth, Location_of_IfcExtrudedAreaSolid, Axis, RefDirection = Item_IfcExtrudedAreaSolid_analysis(Item_IfcExtrudedAreaSolid)
                    # Here we distinguish between IfcElements with LCS Rotated in terms of GCS : IfcArbitraryClosedProfileDef, or parallel with GCS: IfcRectangleProfileDef
                    
                    #Step 1: a) Define Local origin point of Ifc Element in terms of GCS: Local_Orgin_Coord
                    #        b) Define the Vertices of the IfcElement in hypotherical scenario where the LCS is parallel with GCS 
                    print(Local_Origin_Ifcobject)
                    #Here we store the Coordinates of the Local Origin point of the IfcWall in terms of the GCS
                    Local_Orgin_Coord = np.array(Local_Origin_Ifcobject)
                        
                    # we take the SlabPoints variable (IfcPolyline) that was produced before.
                    print(f"SlabPoints", SlabPoints)
                    Point_A = SlabPoints[2].Coordinates
                    Point_A = np.array(Point_A)
                    print("Coordinates of first point B of IFCWall area:", Point_A)
                    Point_B = SlabPoints[1].Coordinates
                    Point_B = np.array(Point_B)
                    print("Coordinates of first point B of IFCWall area:", Point_B)
                    Point_C = SlabPoints[0].Coordinates
                    Point_C = np.array(Point_C)
                    print("Coordinates of first point A of IFCWall area:", Point_C)
                    Point_D = SlabPoints[3].Coordinates
                    Point_D = np.array(Point_D)
                    print("Coordinates of first point B of IFCWall area:", Point_D)
                    #We identify the Z coordinate of the Level in which the IfcSlab is located
                    Coordinates_of_FloorLevel = Location2.Coordinates
                    print(f"Coordinates_of_FloorLevel", Coordinates_of_FloorLevel)
                    Coordinates_of_FloorLevel = np.array(Coordinates_of_FloorLevel)
                    Floor_Z_Coord_slab = Coordinates_of_FloorLevel[2]
                    print(f"Floor_Z_Coord", Floor_Z_Coord_slab)
                    #Then we identify whether or not there is an Height Offset from level, 
                    #for the Extruded area of the Ifcslab from the variable: Location_of_IfcExtrudedAreaSolid
                    Location_of_IfcExtrudedAreaSolid_Coords = Location_of_IfcExtrudedAreaSolid.Coordinates
                    print(f"Location_of_IfcExtrudedAreaSolid_Coords", Location_of_IfcExtrudedAreaSolid_Coords)
                    Location_of_IfcExtrudedAreaSolid_Coords = np.array(Location_of_IfcExtrudedAreaSolid_Coords)
                    Height_Offset = Location_of_IfcExtrudedAreaSolid_Coords[2]
                    print(f"Height_Offset", Height_Offset)
                    #Now we define the position in Z-axis of the extruded area of the IfcSlab that contains 
                    #The A1,B1,C1,D1, Always the A2,B2,C2,D2 will be located above the Extruded Area
                    if Height_Offset != 0:
                        Absolute_Coord_of_ExtrArea = Floor_Z_Coord_slab + Height_Offset
                    else:
                        Absolute_Coord_of_ExtrArea = Floor_Z_Coord_slab
                    #for IfcSlab we can find the absolute coordianres of its vertices in terms of GCS directly 
                    #from the ifcShapeReprsntation2: Foot Print, Curved2d which defines the above points, 
                    #namely the x and y coordinates of the four vertices of IfcSlab.
                    #However we still need info about the z coordinate of each vertex.
                    Vert_Correct_LCS = [(Point_A[0], Point_A[1], Absolute_Coord_of_ExtrArea, 1),
                                    (Point_B[0], Point_B[1], Absolute_Coord_of_ExtrArea, 1),
                                    (Point_C[0], Point_C[1], Absolute_Coord_of_ExtrArea, 1),
                                    (Point_D[0], Point_D[1], Absolute_Coord_of_ExtrArea, 1),
                                    (Point_A[0], Point_A[1], Absolute_Coord_of_ExtrArea - Depth, 1),
                                    (Point_B[0], Point_B[1], Absolute_Coord_of_ExtrArea - Depth, 1),
                                    (Point_C[0], Point_C[1], Absolute_Coord_of_ExtrArea - Depth, 1),
                                    (Point_D[0], Point_D[1], Absolute_Coord_of_ExtrArea -Depth, 1)] #The coordinate z is derived from IfcLocalPlacement of the IfcElement
                    Vert_Correct_LCS = np.array(Vert_Correct_LCS)
                    Vert_Correct_LCS = np.delete(Vert_Correct_LCS, -1, axis=1)
                    print(Vert_Correct_LCS)
                    print("Vert_Hyp_LCS", Vert_Correct_LCS)
                    #The desired face is always the down face of the slab namely the face which is defined by the vertices A2,B2,C2,D2
                    Slab_face =[(Point_A[0], Point_A[1], Absolute_Coord_of_ExtrArea - Depth),
                                    (Point_B[0], Point_B[1], Absolute_Coord_of_ExtrArea - Depth),
                                    (Point_C[0], Point_C[1], Absolute_Coord_of_ExtrArea - Depth),
                                    (Point_D[0], Point_D[1], Absolute_Coord_of_ExtrArea -Depth)] 
                    Slab_face = np.array(Slab_face)
                    print("The desired Slab face is given by these coordinates:", Slab_face)
                    # Creating a list of the matrices, its GUID, and its ID
                    IfcSlab_Matrices.append((Vert_Correct_LCS, Slab_face, GUid, N_Slab, Absolute_Coord_of_ExtrArea ))
                    New_slab_location_object = np.array([ "IfcSlab", GUid, Absolute_Coord_of_ExtrArea ], dtype=object)
                    Slab_Z_Location = np.vstack([Slab_Z_Location, New_slab_location_object])
                    
                else:
                    print("this Slab is not represented by SweptSolid Representation Type, so it is excluded from the analysis")
                if N_Slab == 1:
                    faces = np.array([
                    [0, 1, 5, 4],  # 1st XZ face
                    [1, 2, 6, 5],  # 1st YZ face
                    [3, 0, 1, 2],  # 1st XY face
                    [2, 3, 7, 6],  # 2st XZ face
                    [3, 0, 4, 7],  # 2st YZ face
                    [7, 4, 5, 6]   # 2st XY face
                    ])
                # Separate faces into separate arrays
                separate_faces = [Vert_Correct_LCS[face_indices] for face_indices in faces]
                New_object = np.array([patternIfcSlab, GUid, IfcIdentifier, N_Slab, Absolute_Coord_of_ExtrArea, None, None, None, None, None, None], dtype=object)
                #Iterate through the number of faces 
                for i, face in enumerate(separate_faces):
                    if i==0: #face 1 XZ in Complete matrix
                        New_object[6] = face  # Assign the first value of each array to the corresponding column in new_rowNew_object[i + 5] = face  # Assign the first value of each array to the corresponding column in new_row
                    elif i==1: #Face 2 YZ in Complete matrix
                        New_object[5] = face  # Assign the first value of each array to the corresponding column in new_row
                    elif i==2: #face 3 XY in Complete matrix
                        New_object[7] = face
                    elif i==3: #face 4 XZ in Complete matrix
                        New_object[9] = face
                    elif i==4: #face 5 YZ in Complete matrix
                        New_object[8] = face
                    elif i==5: #face 6 XY in Complete matrix
                        New_object[10] = face
                
                # Append the new row to the array
                Complete_Matrix = np.vstack([Complete_Matrix, New_object])
                print(Complete_Matrix)      
#Here we investigate which slab is located in the least z coordinate. This slab will constitute then the lower celing of the floor, in which every element is located.
#With this way we exclude the Upper slab. This will be used in visualisation part of the algorithm anmely in the part that we create the meshes of the Building Elements and store them in a OBJ file
Z_coords = Slab_Z_Location[1:, 2]
Max_z_coord = max(Z_coords)
for slab in Slab_Z_Location[1:]:
    if slab[2] == Max_z_coord:
        Higher_slab_GUId = slab[1]
# Walls
# File path
ifc_file_path_Wall = r"C:/Python Output/Real data/TXT files/SDK_walls.txt"
with open(ifc_file_path_Wall, "w") as file:
    # Iterate through the list of tuples
    for matrix, identifier, index, z in IfcWall_Matrices:
            file.write(f"Wall {index} ({identifier}) Floor {z}:\n")
            np.savetxt(file, matrix, fmt='%f', delimiter=', ')
            file.write("\n\n")
# File path
ifc_file_path_OpeningElement = r"C:/Python Output/Real data/TXT files/SDK__OpeningElements.txt"
with open(ifc_file_path_OpeningElement, "w") as file:
    # Iterate through the list of tuples
    for matrix, identifier, Tag, index, z in IfcOpeningElement_Matrices:
            file.write(f"Wall {index} ({identifier}) Floor {z}:\n")
            np.savetxt(file, matrix, fmt='%f', delimiter=', ')
            file.write("\n\n")
# Columns
# File path
ifc_file_path_Column = r"C:/Python Output/Real data/TXT files/SDK_Columns.txt"       
with open(ifc_file_path_Column, "w") as file:
    # Iterate through the list of tuples
    for matrix, identifier, index, z in IfcColumn_Matrices:
        file.write(f"Column {index} ({identifier}) Floor {z}:\n")
        np.savetxt(file, matrix, fmt='%f', delimiter=', ')
        file.write("\n\n")
# Beams
# File path
ifc_file_path_Beam = r"C:/Python Output/Real data/TXT files/SDK_Beams.txt"
with open(ifc_file_path_Beam, "w") as file:
    # Iterate through the list of tuples
    for matrix, identifier, index, z in IfcBeam_Matrices:
        file.write(f"Beam {index} ({identifier}) Floor {z}:\n")
        np.savetxt(file, matrix, fmt='%f', delimiter=', ')
        file.write("\n\n")
# Slabs
# File path
ifc_file_path_Slab = r"C:/Python Output/Real data/TXT files/SDK_Slabs.txt"
with open(ifc_file_path_Slab, "w") as file:
    # Iterate through the list of tuples
    for matrix1, matrix2, identifier, index, z in IfcSlab_Matrices:
        file.write(f"Slab {index} ({identifier}) Floor {z}:\n")
        np.savetxt(file, matrix1, fmt='%f', delimiter=', ')
        file.write(f"Desirable face of Slab {index} ({identifier}):\n")
        np.savetxt(file, matrix2, fmt='%f', delimiter=', ')
        file.write("\n\n")
 # File path
def reduce_scale_from_mm_to_m(Complete_Matrix):
    # Iterate through rows
    j = 0
    for row in Complete_Matrix:
        j = j +1
        if j != 1:
            # Iterate through columns 6 to 11
            for i in range(5, 11):
                # Divide values by 1000
                if isinstance(row[i], np.ndarray):
                    row[i] /= 1000
    return Complete_Matrix
Complete_Matrix = reduce_scale_from_mm_to_m(Complete_Matrix)
# File path
ifc_file_path_Complete_Matrix = r"C:/Python Output/Real data/TXT files/SDK_Complete_Matrix.txt"
with open(ifc_file_path_Complete_Matrix, "w") as file:
    # Iterate through the list of tuples
    for IFCtype, guid, identifier, localid, z, face1XZ, face2YZ, face3XY, face4XZ, face5YZ, face6XY in Complete_Matrix[1:]:
            file.write(f" {IFCtype} {guid} {localid} ({identifier}) Floor {z}:\n")
            file.write(f" Face 1 XZ :\n")
            np.savetxt(file, face1XZ, fmt='%f', delimiter=', ')
            file.write(f" Face 2 YZ :\n")
            np.savetxt(file, face2YZ, fmt='%f', delimiter=', ')
            file.write(f" Face 3 XY :\n")
            np.savetxt(file, face3XY, fmt='%f', delimiter=', ')
            file.write(f" Face 4 XZ :\n")
            np.savetxt(file, face4XZ, fmt='%f', delimiter=', ')
            file.write(f" Face 5 YZ :\n")
            np.savetxt(file, face5YZ, fmt='%f', delimiter=', ')
            file.write(f" Face 6 XY :\n")
            np.savetxt(file, face6XY, fmt='%f', delimiter=', ')
            file.write("\n\n")
print(IfcWall_Matrices)    
print("END")        
print(Complete_Matrix)
print(Complete_Matrix[1][4])

#WALL TO WALL INTERSECTION CASES (The same procedure should be followed for Wall to Beam intersections, Wall to Column intersections, but in this research we focus only on Walls)
#FIRST: Create a New Complete matrix that will contain only the wall faces
Wall_Complete_Matrix = np.array([["IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face 1 XZ", "Face 2 YZ", "Face 3 XY", "Face 4 XZ", "Face 5 YZ", "Face 6 XY"]], dtype=object)
def wall_face_features(wall_face):
    x = wall_face[:, 0] #This is slicing operation selects all rows (:) and only the elements in the first column (0)
    y = wall_face[:, 1]
    z = wall_face[:, 2]
    # Calculate differences
    dx = np.max(x) - np.min(x)
    dy = np.max(y) - np.min(y)
    dz = np.max(z) - np.min(z)
    return dx, dy, dz, x, y, z
        
for row in Complete_Matrix[1:]:
 if re.search(patternIfcWall, str(row[0])): 
    New_Wall_object = np.array([ row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10]], dtype=object)
    Wall_Complete_Matrix = np.vstack([Wall_Complete_Matrix, New_Wall_object])
print(Wall_Complete_Matrix) #Then we should include them back in the complete matrix
#Array of YZ wall faces
wall_faces_YZ = np.array([["IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face Type in GCS", "Number of Face", "Face YZ", "Average Location", "Other Dimension", "Max", "Min", "Max_z", "Min_z"]], dtype=object)
for row_YZ in Wall_Complete_Matrix[1:]:
    k = 0
    for k in range(2):
        k = k +1
        if k == 1:
            dx, dy, dz, x, y, z = wall_face_features(row_YZ[6])
            if -0.001 <= dy <= 0.001:
                Face_Type_GCS = "Plane XZ" #This face name is refered to the GCS of the project in order to communicate efficiently with the corresponding segment
                AverageLocationY = y[0] #all points in the face have the same y value so it is not neccessary to calculate the averge
                New_Wall_YZ_face = np.array([ row_YZ[0], row_YZ[1], row_YZ[2], row_YZ[3], row_YZ[4], Face_Type_GCS, "Face 1 YZ", row_YZ[6], y[0], dx, max(x), min(x), max(z), min(z)], dtype=object)
                wall_faces_YZ = np.vstack([wall_faces_YZ, New_Wall_YZ_face])
            elif -0.001 <= dx <= 0.001:
                Face_Type_GCS = "Plane YZ"
                AverageLocationX = x[0] #all points in the face have the same x value so it is not necessary to calculate the averge
                New_Wall_YZ_face = np.array([ row_YZ[0], row_YZ[1], row_YZ[2], row_YZ[3], row_YZ[4], Face_Type_GCS, "Face 1 YZ", row_YZ[6], x[0], dy, max(y), min(y), max(z), min(z)], dtype=object)
                wall_faces_YZ = np.vstack([wall_faces_YZ, New_Wall_YZ_face])
        elif k == 2:
            dx, dy, dz, x, y, z = wall_face_features(row_YZ[9])
            if -0.001 <= dy <= 0.001:
                Face_Type_GCS = "Plane XZ" #This face name is refered to the GCS of the project in order to communicate efficiently with the corresponding segment
                AverageLocationY = y[0] #all points in the face have the same y value so it is not neccessary to calculate the averge
                New_Wall_YZ_face = np.array([ row_YZ[0], row_YZ[1], row_YZ[2], row_YZ[3], row_YZ[4], Face_Type_GCS, "Face 2 YZ", row_YZ[9], y[0], dx, max(x), min(x), max(z), min(z)], dtype=object)
                wall_faces_YZ = np.vstack([wall_faces_YZ, New_Wall_YZ_face])
            elif -0.001 <= dx <= 0.001:
                Face_Type_GCS = "Plane YZ"
                AverageLocationX = x[0] #all points in the face have the same x value so it is not necessary to calculate the averge
                New_Wall_YZ_face = np.array([ row_YZ[0], row_YZ[1], row_YZ[2], row_YZ[3], row_YZ[4], Face_Type_GCS, "Face 2 YZ", row_YZ[9], x[0], dy, max(y), min(y), max(z), min(z)], dtype=object)
                wall_faces_YZ = np.vstack([wall_faces_YZ, New_Wall_YZ_face])
print(wall_faces_YZ)
#Now we have all the YZ faces that we want to check in an array
#Let's put also the XZ faces that we want to check in an array
#Array of XZ wall faces. Be careful the XZ, and YZ notation refers to the LCS of the corresponding walls.
wall_faces_XZ = np.array([["IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face Type in GCS", "Number of Face", "Face XZ", "Average Location", "Other Dimension", "Max", "Min", "Max_z", "Min_z"]], dtype=object)
for row_XZ in Wall_Complete_Matrix[1:]:
    k = 0
    for k in range(2):
        k = k +1
        if k == 1:
            dx, dy, dz, x, y, z = wall_face_features(row_XZ[5])
            if -0.001 <= dy <= 0.001:
                Face_Type_GCS = "Plane XZ" #This face name is refered to the GCS of the project in order to communicate efficiently with the corresponding segment
                AverageLocationY = y[0] #all points in the face have the same y value so it is not neccessary to calculate the averge
                New_Wall_XZ_face = np.array([ row_XZ[0], row_XZ[1], row_XZ[2], row_XZ[3], row_XZ[4], Face_Type_GCS, "Face 1 XZ", row_XZ[5], y[0], dx, max(x), min(x), max(z), min(z)], dtype=object)
                wall_faces_XZ = np.vstack([wall_faces_XZ, New_Wall_XZ_face])
            elif -0.001 <= dx <= 0.001:
                Face_Type_GCS = "Plane YZ"
                AverageLocationX = x[0] #all points in the face have the same x value so it is not necessary to calculate the averge
                New_Wall_XZ_face = np.array([ row_XZ[0], row_XZ[1], row_XZ[2], row_XZ[3], row_XZ[4], Face_Type_GCS, "Face 1 XZ", row_XZ[5], x[0], dy, max(y), min(y), max(z), min(z)], dtype=object)
                wall_faces_XZ = np.vstack([wall_faces_XZ, New_Wall_XZ_face])
        elif k == 2:
            dx, dy, dz, x, y, z = wall_face_features(row_XZ[8])
            if -0.001 <= dy <= 0.001:
                Face_Type_GCS = "Plane XZ" #This face name is refered to the GCS of the project in order to communicate efficiently with the corresponding segment
                AverageLocationY = y[0] #all points in the face have the same y value so it is not neccessary to calculate the averge
                New_Wall_XZ_face = np.array([ row_XZ[0], row_XZ[1], row_XZ[2], row_XZ[3], row_XZ[4], Face_Type_GCS, "Face 2 XZ", row_XZ[8], y[0], dx, max(x), min(x), max(z), min(z)], dtype=object)
                wall_faces_XZ = np.vstack([wall_faces_XZ, New_Wall_XZ_face])
            elif -0.001 <= dx <= 0.001:
                Face_Type_GCS = "Plane YZ"
                AverageLocationX = x[0] #all points in the face have the same x value so it is not necessary to calculate the averge
                New_Wall_XZ_face = np.array([ row_XZ[0], row_XZ[1], row_XZ[2], row_XZ[3], row_XZ[4], Face_Type_GCS, "Face 2 XZ", row_XZ[8],  x[0], dy, max(y), min(y), max(z), min(z)], dtype=object)
                wall_faces_XZ = np.vstack([wall_faces_XZ, New_Wall_XZ_face])
print(wall_faces_XZ)
#case
# Process of modifying the Vertcies of a Face_XZ of a wall
# First we create an array that we will store all the adjusted faces, and their attributes.
Adjusted_wall_faces_XZ = np.array([["IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face Type in GCS", "Number of Face", "Face XZ", "Average Location", "Other Dimension", "Max", "Min", "Max_z", "Min_z"]], dtype=object)

array_XZ = None
array_YZ = None
for wall_face_YZ in wall_faces_YZ[1:]:
    for wall_face_XZ in wall_faces_XZ[1:]:
        #if wall_face_XZ[1] == "1mNm2D_bLEcw7ElQBoHPDF" and wall_face_YZ[1] == "1mNm2D_bLEcw7ElQBoHPD0":
        #if wall_face_XZ[1] == "1mNm2D_bLEcw7ElQBoHOn9" and wall_face_YZ[1] == "1mNm2D_bLEcw7ElQBoHOn8":
        #if wall_face_XZ[1] == "1mNm2D_bLEcw7ElQBoHPJa" and wall_face_YZ[1] == "1mNm2D_bLEcw7ElQBoHPJb":
        #if wall_face_XZ[1] == "1mNm2D_bLEcw7ElQBoHPJa" and wall_face_YZ[1] == "1mNm2D_bLEcw7ElQBoHPD1":
        #if wall_face_XZ[1] == "1mNm2D_bLEcw7ElQBoHPD_" and wall_face_YZ[1] == "1mNm2D_bLEcw7ElQBoHPCH":
            #Condition 1:
            condit_1 = wall_face_YZ[8] - wall_face_XZ[8]
            #Condition 3.1: Max
            condit_3_1 = wall_face_YZ[10] - wall_face_XZ[10]
            #Condition 3.2: Min
            condit_3_2 = wall_face_YZ[11] - wall_face_XZ[11]
            if (-0.001 <= condit_1 <= 0.001 and wall_face_YZ[5] == wall_face_XZ[5] and (-0.001 <= condit_3_1 <= 0.001 or -0.001 <= condit_3_2 <= 0.001)):
                #Then the conditions are satisfied so we can modify the vertices of the XZ face.
                if -0.001 <= condit_3_1 <= 0.001: #the two faces have the same Max coordinate
                    array_YZ = wall_face_YZ[7]
                    array_XZ = wall_face_XZ[7]
                    if wall_face_YZ[5] == "Plane XZ": #then the two faces belong to the XZ plane in terms of GCS so they both have common coordinate y in the second column of their face array.
                        #Becasue we have the same max coordinate between those two face we cannot compare them with this. So we compare them with their Min x coordinates.
                        #Becasue we are in the plane XZ in terms of GCS we capture the MIN x coordinate of face YZ. 
                        #This Min x coordinate of face YZ will substitute the Min x Coordinate of the Face XZ.
                        #Min and Max Coordinate of Face YZ
                        min_YZ = wall_face_YZ[11] # this will substitute the min coordinate of face_XZ
                        max_YZ = wall_face_YZ[10]
                        #we check the Array XZ
                        x = array_XZ[:, 0]
                        for i in range(len(x)):
                            x_coord = array_XZ[i][0] #we want to check the Array XZ, and the coordinate y so the index 1. This variable might be the min or the max of the y coordinate of Face_XZ
                            #the process continues like this: if the y_coord != max_YZ then this coordinate represents a mean and we substitute this mean with the min_YZ 
                            if round(x_coord, 3) == round(max_YZ, 3):
                                array_XZ[i][0] = min_YZ
                        x_new = array_XZ[:, 0]
                        dx = np.max(x_new) - np.min(x_new)
                        new_adjusted_XZ_object = np.array([ wall_face_XZ[0], wall_face_XZ[1], wall_face_XZ[2], wall_face_XZ[3], wall_face_XZ[4], wall_face_XZ[5], wall_face_XZ[6], array_XZ, wall_face_XZ[8], dx, max(x_new), wall_face_XZ[11], wall_face_XZ[12], wall_face_XZ[13]], dtype=object)
                        Adjusted_wall_faces_XZ = np.vstack([Adjusted_wall_faces_XZ, new_adjusted_XZ_object])
                    if wall_face_YZ[5] == "Plane YZ":  #then the two faces belong to the YZ plane in terms of GCS so they both have common coordinate x in the first column of their face array.
                        #Becasue we have the same max coordinate between those two face we cannot compare them with this. So, again, we compare them with their Min y coordinates.
                        #Becasue we are in the plane XZ in terms of GCS we capture the MIN y coordinate of face YZ. 
                        #This Min y coordinate of face YZ will substitute the Min y Coordinate of the Face XZ.
                        #Min and Max Coordinate of Face YZ
                        min_YZ = wall_face_YZ[11] # this will substitute the min coordinate of face_XZ
                        max_YZ = wall_face_YZ[10]
                        #we check the Array XZ
                        y = array_XZ[:, 1] 
                        for i in range(len(y)):
                            y_coord = array_XZ[i][1] #we want to check the Array XZ, and the coordinate y so the index 1. This variable might be the min or the max of the y coordinate of Face_XZ
                            #the process continues like this: if the y_coord != max_YZ then this coordinate represents a mean and we substitute this mean with the min_YZ 
                            if round(y_coord, 3) == round(max_YZ, 3):
                                array_XZ[i][1] = min_YZ
                        y_new = array_XZ[:, 0]
                        dy = np.max(y_new) - np.min(y_new)
                        new_adjusted_XZ_object = np.array([ wall_face_XZ[0], wall_face_XZ[1], wall_face_XZ[2], wall_face_XZ[3], wall_face_XZ[4], wall_face_XZ[5], wall_face_XZ[6], array_XZ, wall_face_XZ[8], dy, max(y_new), wall_face_XZ[11], wall_face_XZ[12], wall_face_XZ[13]], dtype=object)
                        Adjusted_wall_faces_XZ = np.vstack([Adjusted_wall_faces_XZ, new_adjusted_XZ_object])
                elif -0.001 <= condit_3_2 <= 0.001: #the two faces have the same Min coordinate
                        array_YZ = wall_face_YZ[7]
                        array_XZ = wall_face_XZ[7]
                        if wall_face_YZ[5] == "Plane XZ": #then the two faces belong to the XZ plane in terms of GCS so they both have common coordinate y in the second column of their face array.
                            # Becasue we have the same Min coordinate between those two face we cannot compare them with this. So we compare them with their Max x coordinates.
                            # Becasue we are in the plane XZ in terms of GCS we capture the Max x coordinate of face YZ. 
                            # This Max x coordinate of face YZ will substitute the Max x Coordinate of the Face XZ.
                            # Min and Max Coordinate of Face YZ
                            min_YZ = wall_face_YZ[11] 
                            max_YZ = wall_face_YZ[10] # this will substitute the max coordinate of face_XZ
                            #we check the Array XZ
                            x = array_XZ[:, 0]
                            for i in range(len(x)):
                                x_coord = array_XZ[i][0] #we want to check the Array XZ, and the coordinate x so the index 0. This variable might be the min or the max of the x coordinate of Face_XZ
                                # the process continues like this: if the y_coord != min_YZ then this coordinate represents a max and we substitute this max with the max_YZ 
                                if round(x_coord, 3) == round(min_YZ, 3):
                                    array_XZ[i][0] = max_YZ
                            x_new = array_XZ[:, 0]
                            dx = np.max(x_new) - np.min(x_new)
                            new_adjusted_XZ_object = np.array([ wall_face_XZ[0], wall_face_XZ[1], wall_face_XZ[2], wall_face_XZ[3], wall_face_XZ[4], wall_face_XZ[5], wall_face_XZ[6], array_XZ, wall_face_XZ[8], dx, max(x_new), wall_face_XZ[11], wall_face_XZ[12], wall_face_XZ[13]], dtype=object)
                            Adjusted_wall_faces_XZ = np.vstack([Adjusted_wall_faces_XZ, new_adjusted_XZ_object])
                        if wall_face_YZ[5] == "Plane YZ":  #then the two faces belong to the YZ plane in terms of GCS so they both have common coordinate x in the first column of their face array.
                            #Becasue we have the same Min coordinate between those two face we cannot compare them with this. So, again, we compare them with their Max y coordinates.
                            #Becasue we are in the plane XZ in terms of GCS we capture the MIN y coordinate of face YZ. 
                            #This Max y coordinate of face YZ will substitute the Max y Coordinate of the Face XZ.
                            #Min and Max Coordinate of Face YZ
                            min_YZ = wall_face_YZ[11] 
                            max_YZ = wall_face_YZ[10] # this will substitute the min coordinate of face_XZ
                            #we check the Array XZ
                            y = array_XZ[:, 1] 
                            for i in range(len(y)):
                                y_coord = array_XZ[i][1] #we want to check the Array XZ, and the coordinate x so the index 0. This variable might be the min or the max of the x coordinate of Face_XZ
                                # the process continues like this: if the y_coord != min_YZ then this coordinate represents a max and we substitute this max with the max_YZ 
                                if round(y_coord, 3) == round(min_YZ, 3):
                                    array_XZ[i][1] = max_YZ
                            y_new = array_XZ[:, 0]
                            dy = np.max(y_new) - np.min(y_new)
                            new_adjusted_XZ_object = np.array([ wall_face_XZ[0], wall_face_XZ[1], wall_face_XZ[2], wall_face_XZ[3], wall_face_XZ[4], wall_face_XZ[5], wall_face_XZ[6], array_XZ, wall_face_XZ[8], dy, max(y_new), wall_face_XZ[11], wall_face_XZ[12], wall_face_XZ[13]], dtype=object)
                            Adjusted_wall_faces_XZ = np.vstack([Adjusted_wall_faces_XZ, new_adjusted_XZ_object])
print(Adjusted_wall_faces_XZ)
#We have found the new, adjusted wall faces. Now we have to substitute them in the complete matrix
#We Iterate over the variable Adjusted_wall_faces_XZ, and then we check if it equals with a wall in the Complete_matrrix array. if it equals then we suhbstitute its correspondig face

for adjusted_face in Adjusted_wall_faces_XZ[1:]:
    for index, Primary_face in enumerate(Complete_Matrix[1:], start=1):
        if adjusted_face[1] == Primary_face[1]: #If the GUIDs between the two walls from the two faces are the same then we check which face do we have
            if adjusted_face[6] == "Face 1 XZ":
                Complete_Matrix[index][5] = adjusted_face[7] # We substitute the previous face with the adjusted face
            elif adjusted_face[6] == "Face 2 XZ":
                Complete_Matrix[index][8] = adjusted_face[7] # We substitute the previous face with the adjusted face
print(Complete_Matrix)


#COMPARISON PROCEDURE
#STEP 1: Create arrays Face Xz, Face YZ, Face Xy to store, for each face in the model, all the info necessary for the comparison 
individual_triangle = np.array([
    [0, 1, 2 ],  # Tr1 
    [2, 3, 0 ],  # Tr2 
    ])
# Create a list to store mesh objects
# Iterating over the faces
normals = []
Complete_Normal_Matrix = np.array([["IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face 1 XZ", "Face 2 YZ", "Face 3 XY", "Face 4 XZ", "Face 5 YZ", "Face 6 XY","Normal F1","Normal F2", "Nomral F3", "Nomral F4","Nomral F5", "Nomral F6" ]], dtype=object)
FacesXZ = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin" ]], dtype=object)
FacesYZ = np.array([["Face Type in GCS", "Face Type in LCS","IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin" ]], dtype=object)
FacesXY = np.array([["Face Type in GCS", "Face Type in LCS","IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Ymax", "Ymin" ]], dtype=object)
FacesXZ_OE = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin" ]], dtype=object)
FacesYZ_OE = np.array([["Face Type in GCS", "Face Type in LCS","IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin" ]], dtype=object)
FacesXY_OE = np.array([["Face Type in GCS", "Face Type in LCS","IfcType", "GUID", "ID", "Local ID", "Floor Level (z)", "Face Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Ymax", "Ymin" ]], dtype=object)
AverageLocations = [] #we create this array in order to be able to verify the results in the below matrices

def Face_Type_In_terms_of_LCS(face_i):
    if face_i == 1:
        Face_Type_LCS = "Face 1 XZ"
    elif face_i == 2:
        Face_Type_LCS = "Face 1 YZ"
    elif face_i == 3:
        Face_Type_LCS = "Face 1 XY"
    elif face_i == 4:
        Face_Type_LCS = "Face 2 XZ"
    elif face_i == 5:
        Face_Type_LCS = "Face 2 YZ"
    elif face_i == 6:
        Face_Type_LCS = "Face 2 XY"
    return Face_Type_LCS
#Calculation of Complete Matrix with normals, and Face XZ, YZ, XY arrays in terms of GCS
def Calculations_normals_and_faces_in_terms_of_GCS(row, FacesXZ, FacesYZ, FacesXY):
    New_object = np.array([row[0], row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], None, None, None, None, None, None], dtype=object)
    j = 10
    face_i=0
    for face in row[5:]:
        j = j +1
        face_i= face_i + 1
        if face_i==3: # This index refers to Face 3 XY, and forces the desirable face normal to always look down, which is the preferable direction for this face
            v1 = face[0]  #-1 specifies the real index in the array as python start count from 0 not 1
            v2 = face[1] 
            v3 = face[2] 
            face_v1 = v2 - v1
            face_v2 = v3 - v1
            normal = np.cross(face_v2, face_v1)
            normal /= np.linalg.norm(normal)  # Normalize
            normals.append(normal)
            New_object[j]= normal
        else:
            # Calculate normals for each face
            v1 = face[0]  #-1 specifies the real index in the array as python start count from 0 not 1
            v2 = face[1] 
            v3 = face[2] 
            face_v1 = v2 - v1
            face_v2 = v3 - v1
            normal = np.cross(face_v1, face_v2)
            normal /= np.linalg.norm(normal)  # Normalize
            normals.append(normal)
            New_object[j]= normal
        # Extract x, y, and z coordinates
            
        x = face[:, 0] #This is slicing operation selects all rows (:) and only the elements in the first column (0)
        y = face[:, 1]
        z = face[:, 2]

        # Calculate differences
        dx = np.max(x) - np.min(x)
        dy = np.max(y) - np.min(y)
        dz = np.max(z) - np.min(z)
        if -0.001 <= dy <= 0.001:
            Face_Type_GCS = "Face XZ" #This face name is refered to the GCS of the project in order to communicate efficiently with the corresponding segment
            Face_Type_LCS = Face_Type_In_terms_of_LCS(face_i) # This face name is refered to the LCS of each IfcObject in which the face belongs, this is to know in which face we refer to every time
            AverageLocationY = y[0] #all points in the face have the same y value so it is not neccessary to calculate the avergae
            AverageLocationX = (min(x) + max(x))/2
            AverageLocationZ = (min(z) + max(z))/2
            New_face_object = np.array([Face_Type_GCS, Face_Type_LCS, row[0], row[1], row[2], row[3], row[4], y[0], dx, dz, AverageLocationX, AverageLocationZ, normal, face, max(x), min(x), max(z), min(z) ], dtype=object) #the Average Locations will be used later on Clashing Faces with Segments
            FacesXZ = np.vstack([FacesXZ, New_face_object])
            AverageLocations.append(AverageLocationY)
        elif -0.001 <= dx <= 0.001:
            Face_Type_GCS = "Face YZ"
            Face_Type_LCS = Face_Type_In_terms_of_LCS(face_i)
            AverageLocationX = x[0] #all points in the face have the same x value so it is not necessary to calculate the avergae
            AverageLocationY = (min(y) + max(y))/2
            AverageLocationZ = (min(z) + max(z))/2
            New_face_object = np.array([Face_Type_GCS, Face_Type_LCS, row[0], row[1], row[2], row[3], row[4], x[0], dy, dz, AverageLocationY, AverageLocationZ, normal, face, max(y), min(y), max(z), min(z)], dtype=object)
            FacesYZ = np.vstack([FacesYZ, New_face_object])
            AverageLocations.append(AverageLocationX)
        elif -0.001 <= dz <= 0.001:
            Face_Type_GCS = "Face XY"
            Face_Type_LCS = Face_Type_In_terms_of_LCS(face_i)
            AverageLocationZ = z[0] #all points in the face have the same z value so it is not necessary to calculate the avergae
            AverageLocationX = (min(x) + max(x))/2
            AverageLocationY = (min(y) + max(y))/2
            New_face_object = np.array([Face_Type_GCS, Face_Type_LCS, row[0], row[1], row[2], row[3], row[4], z[0], dx, dy, AverageLocationX, AverageLocationY, normal, face, max(x), min(x), max(y), min(y)], dtype=object)
            FacesXY = np.vstack([FacesXY, New_face_object])
            AverageLocations.append(AverageLocationZ)
    return FacesXZ, FacesYZ, FacesXY, AverageLocations, normals, New_object
#By using the above function we manage to produce 1) the arrays FacesXZ, FacesYZ, FacesXY, that contain the faces of Every Building element, and which are used in the Comparison procedure
# and 2) the arrays FacesXZ_OE, FacesYZ_OE, FacesXY_OE, that contain the faces only of the Opening Elements and will be used in the second compariosn procedure which refers only to opening elements
for row in Complete_Matrix[1:]:
    if re.search(patternIfcOpeningElement, str(row[0])):
        FacesXZ_OE, FacesYZ_OE, FacesXY_OE, AverageLocations, normals, New_object = Calculations_normals_and_faces_in_terms_of_GCS(row, FacesXZ_OE, FacesYZ_OE, FacesXY_OE)
    else:
        FacesXZ, FacesYZ, FacesXY, AverageLocations, normals, New_object = Calculations_normals_and_faces_in_terms_of_GCS(row, FacesXZ, FacesYZ, FacesXY)
        Complete_Normal_Matrix = np.vstack([Complete_Normal_Matrix, New_object])
AverageLocationsNP = np.array(AverageLocations)
print(normals)
print(AverageLocationsNP)

#STEP 2: Create arrays Segment XZ, Segment YZ, Segment XY, to store, for each segment, all the info necessary for the comparison
#Analyse .ply files in loction: Segments_path = r"C:/Python Output/Testdata/PLY files/GT_Topoplan_psich/Column1_only/Segments/PLY_GT_Topoplan_2_psich_Cube_selected_Cube_Column1_Only_Translated_Rotated.ply"
# Directory where your .ply files are located
input_Segement_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_final/Segments/'
Segement_XZ= np.array([["Segment ID", "Segment Name", "Segment Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Xmax", "Xmin", "Zmax", "Zmin", "Colour" ]], dtype=object)
Segement_YZ= np.array([["Segment ID", "Segment Name", "Segment Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Ymax", "Ymin", "Zmax", "Zmin", "Colour" ]], dtype=object)
Segement_XY= np.array([["Segment ID", "Segment Name", "Segment Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Xmax", "Xmin", "Ymax", "Ymin", "Colour" ]], dtype=object)
output_Segment_directory = 'C:/Python Output/Real Data/PLY files/SDK_real_dataset_final/Segments_New/'
# Iterate over files in the directory
Seg_XZ_i = 0
Seg_YZ_i = 0
Seg_XY_i = 0
def Rewrite_segement_ply_files(Segment_Name, Seg_i, output_directory):
    # Save the processed point cloud to a new file
    output_filename = f"{Segment_Name}_{Seg_i}.ply"
    output_filepath = os.path.join(output_directory, output_filename)
    os.makedirs(output_directory, exist_ok=True)
    o3d.io.write_point_cloud(output_filepath, ply)

for filename in os.listdir(input_Segement_directory):
    if filename.endswith('.ply'):
        # Load .ply file
        filepath = os.path.join(input_Segement_directory, filename)
        ply = o3d.io.read_point_cloud(filepath)
        # Get points as numpy array
        points_array = np.asarray(ply.points)
        colour_array = np.asarray(ply.colors)
        Seg_x = points_array[:, 0] #This is slicing operation selects all rows (:) and only the elements in the first column (0)
        Seg_y = points_array[:, 1]
        Seg_z = points_array[:, 2]

        # Calculate differences
        Seg_dx = np.max(Seg_x) - np.min(Seg_x)
        Seg_dy = np.max(Seg_y) - np.min(Seg_y)
        Seg_dz = np.max(Seg_z) - np.min(Seg_z)
        # Visualize point cloud
        if Seg_dy < Seg_dx and Seg_dy < Seg_dz:
            Seg_XZ_i = Seg_XZ_i +1
            Segment_Name = "Segment XZ" #This face name is refered to the GCS of the project in order to communicate with the corresponding segment
            average_value = np.mean(Seg_y)
            AverageLocationSegmentX = (min(Seg_x) + max(Seg_x))/2
            AverageLocationSegmentZ = (min(Seg_z) + max(Seg_z))/2
            New_Segment_object = np.array([Seg_XZ_i, Segment_Name, average_value, Seg_dx, Seg_dz, AverageLocationSegmentX, AverageLocationSegmentZ, max(Seg_x), min(Seg_x), max(Seg_z), min(Seg_z), colour_array[0]], dtype=object)
            Segement_XZ = np.vstack([Segement_XZ, New_Segment_object])
            Rewrite_segement_ply_files(Segment_Name, Seg_XZ_i, output_Segment_directory)
        elif Seg_dx < Seg_dy and Seg_dx < Seg_dz:
            Seg_YZ_i = Seg_YZ_i +1
            Segment_Name = "Segment YZ"
            average_value = np.mean(Seg_x)
            AverageLocationSegmentY = (min(Seg_y) + max(Seg_y))/2
            AverageLocationSegmentZ = (min(Seg_z) + max(Seg_z))/2
            New_Segment_object = np.array([Seg_YZ_i, Segment_Name, average_value, Seg_dy, Seg_dz, AverageLocationSegmentY, AverageLocationSegmentZ, max(Seg_y), min(Seg_y), max(Seg_z), min(Seg_z), colour_array[0]], dtype=object)
            Segement_YZ = np.vstack([Segement_YZ, New_Segment_object])
            Rewrite_segement_ply_files(Segment_Name, Seg_YZ_i, output_Segment_directory)
        elif Seg_dz < Seg_dx and Seg_dz < Seg_dy:
            Seg_XY_i = Seg_XY_i +1
            Segment_Name = "Segment XY"
            average_value = np.mean(Seg_z)
            AverageLocationSegmentX = (min(Seg_x) + max(Seg_x))/2
            AverageLocationSegmentY = (min(Seg_y) + max(Seg_y))/2
            New_Segment_object = np.array([Seg_XY_i, Segment_Name, average_value, Seg_dx, Seg_dy, AverageLocationSegmentX, AverageLocationSegmentY, max(Seg_x), min(Seg_x), max(Seg_y), min(Seg_y), colour_array[0]], dtype=object)
            Segement_XY = np.vstack([Segement_XY, New_Segment_object])
            Rewrite_segement_ply_files(Segment_Name, Seg_XY_i, output_Segment_directory)
            
#Check for conflicts or clashes between ifcWall, IfcColumn, IfcBeam, and IfcSlab faces in XZ plane
def Delete_clashing_Faces_in_XZandYZ(cleaned_face_array, Excluded_face_array, Original_Face_Array): # the correct face_array is the array which contains the faces that are correct and coloured with green  8kkkkkkkk9ii, 
    indices_to_delete = []
    Cleaned_Original_Face_Array = Original_Face_Array.copy()
    i = 0
    for Primary_row in Original_Face_Array[1:]:
        j = 0
        i = i + 1
        for Secondary_row in Original_Face_Array[1:]:
            j = j + 1
            #if Primary_row[3] == "05hIZEszH5GBM5fSI338vb" and Secondary_row[3] == "05hIZEszH5GBM5fSI338vU":
            #if Primary_row[3] == "1mNm2D_bLEcw7ElQBoHPD_" and Secondary_row[3] == "1mNm2D_bLEcw7ElQBoHPCH":
            #if Primary_row[3] == "1mNm2D_bLEcw7ElQBoHPCH" and Secondary_row[3] == "1mNm2D_bLEcw7ElQBoHPD_":
            #if Primary_row[3] == "05hIZEszH5GBM5fSI338va" and Secondary_row[3] == "05hIZEszH5GBM5fSI338vU":
            Location_difference = abs(abs(Primary_row[7]) - abs(Secondary_row[7]))
            if -0.005 <= Location_difference <= 0.005 and Secondary_row[15] <= Primary_row[10]<= Secondary_row[14]: # check if the faces are belong in the same Location y, to account for mistakes in IFC model, we give a tolerance in diffrence between y's in amount of 0.005 m, namely 5mm.secondly check if the two elements are really attatched or are just located in the same Axis 1 average location but in different Axis 2 location, namely they are not attatched.
                if Primary_row[2] == "IfcWall":
                    if Secondary_row[2] == "IfcColumn": #Clash between wall face and Column face --> we delete the Wall face (Primary row)
                        #We calculate the area of the two faces. Then The face with the smallest area is deleted
                        Wall_area = Primary_row[8]*Primary_row[9]
                        Column_area = Secondary_row[8]*Secondary_row[9]
                        if Column_area > Wall_area: #then the Wall is folded inside the Column so the wall face is excluded
                            cleaned_face_array = np.vstack([cleaned_face_array, Secondary_row])
                            Excluded_face_array = np.vstack([Excluded_face_array, Primary_row])
                            # Delete the row with index i
                            indices_to_delete.append(i)
                        elif Column_area < Wall_area: #then the Column is folded inside the wall so the Column face is excluded
                            cleaned_face_array = np.vstack([cleaned_face_array, Primary_row])
                            Excluded_face_array = np.vstack([Excluded_face_array, Secondary_row])
                            # Delete the row with index j
                            indices_to_delete.append(j)
                    if i != j:
                        if Secondary_row[2] == "IfcWall": #Clash between wall face and Wall face --> we delete the Wall face that is folded to the big one (Primary row)
                            Wall_area_Primary = Primary_row[8]*Primary_row[9]
                            Wall_area_Secondary = Secondary_row[8]*Secondary_row[9]
                            if Wall_area_Secondary > Wall_area_Primary: #then the Primary Wall is folded inside Secondary wall so we exclude the primary wall excluded
                                # Check if the number is in the array
                                if i in indices_to_delete:
                                    print(f"{i} is in the array")
                                else:
                                    print(f"{i} is not in the array")
                                    cleaned_face_array = np.vstack([cleaned_face_array, Secondary_row])
                                    Excluded_face_array = np.vstack([Excluded_face_array, Primary_row])
                                    # Delete the row with index i
                                    indices_to_delete.append(i)
                            elif Wall_area_Secondary < Wall_area_Primary: #then the Secondary Wall is folded inside Primary wall so we exclude the Secondary wall excluded
                                # Check if the number is in the array
                                if j in indices_to_delete:
                                    print(f"{j} is in the array")
                                else:
                                    print(f"{j} is not in the array")
                                    cleaned_face_array = np.vstack([cleaned_face_array, Primary_row])
                                    Excluded_face_array = np.vstack([Excluded_face_array, Secondary_row])
                                    # Delete the row with index j
                                    indices_to_delete.append(j)
                elif Primary_row[2] == "IfcColumn":
                    if Secondary_row[2] == "IfcSlab": #Clash between Column face and Slab face --> we delete the Column face (Primary row)
                        cleaned_face_array = np.vstack([cleaned_face_array, Secondary_row])
                        Excluded_face_array = np.vstack([Excluded_face_array, Primary_row])
                        # Delete the row with index i
                        indices_to_delete.append(i)
                    elif Secondary_row[2] == "IfcBeam": #Clash between  Column face and Beam face --> we delete the Beam face (Primary row)
                        cleaned_face_array = np.vstack([cleaned_face_array, Primary_row])
                        Excluded_face_array = np.vstack([Excluded_face_array, Secondary_row])
                        # Delete the row with index j
                        indices_to_delete.append(j)
                elif Primary_row[2] == "IfcBeam":
                    if Secondary_row[2] == "IfcSlab": #Clash between Beam face and Slab face --> we delete the Beam face (Primary row)
                        cleaned_face_array = np.vstack([cleaned_face_array, Secondary_row])
                        Excluded_face_array = np.vstack([Excluded_face_array, Primary_row])
                        # Delete the row with index i
                        indices_to_delete.append(i)
                    
    Cleaned_Original_Face_Array = np.delete(Cleaned_Original_Face_Array, indices_to_delete, axis=0)
    return Cleaned_Original_Face_Array, Excluded_face_array, indices_to_delete
#Check for conflicts or clashes between ifcWall, IfcColumn, IfcBeam, and IfcSlab faces in XZ plane
def Delete_clashing_Faces_in_XYplane(cleaned_face_array, Excluded_face_array, Original_Face_Array):
    indices_to_delete = []
    Cleaned_Original_Face_Array = Original_Face_Array.copy()
    i = 0
    for Primary_row in Original_Face_Array[1:]:
        j = 0
        i = i + 1
        for Secondary_row in Original_Face_Array[1:]:
            j = j + 1
            Location_difference = abs(abs(Primary_row[7]) - abs(Secondary_row[7]))
            if -0.005 <= Location_difference <= 0.005: # check if the faces are belong in the same Location y, to account for mistakes in IFC model, we give a tolerance in diffrence between y's in amount of 0.005 m, namely 5mm
                if Primary_row[2] == "IfcWall":
                    if Secondary_row[2] == "IfcSlab": #Clash between wall face and Slab face --> we delete the Wall face (Primary row)
                        cleaned_face_array = np.vstack([cleaned_face_array, Secondary_row])
                        Excluded_face_array = np.vstack([Excluded_face_array, Primary_row])
                        # Delete the row with index i
                        indices_to_delete.append(i)
                elif Primary_row[2] == "IfcColumn":
                    if Secondary_row[2] == "IfcSlab": #Clash between Column face and Slab face --> we delete the Column face (Primary row)
                        cleaned_face_array = np.vstack([cleaned_face_array, Secondary_row])
                        Excluded_face_array = np.vstack([Excluded_face_array, Primary_row])
                        # Delete the row with index i
                        indices_to_delete.append(i)
                elif Primary_row[2] == "IfcBeam": # This part also accounts for the case in which the down face of a beam is in the same location as the down face of slab, then the beam face is considered correct and is not included in the analysis.
                    if Secondary_row[2] == "IfcSlab": #and check the case for which Slabs down XY face is in the same position with Beams down XY face #Clash between Beam face and Slab face --> we delete the Beam face (Primary row)
                        cleaned_face_array = np.vstack([cleaned_face_array, Secondary_row])
                        Excluded_face_array = np.vstack([Excluded_face_array, Primary_row])
                        # Delete the row with index i
                        indices_to_delete.append(i)
    Cleaned_Original_Face_Array = np.delete(Cleaned_Original_Face_Array, indices_to_delete, axis=0)
    return Cleaned_Original_Face_Array, Excluded_face_array, indices_to_delete

# Here we check only if there is a clash between the two faces of different IfcObjects, and we delete one of these faces.
Cleaned_FaceXZ_inp= np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin"]], dtype=object)
Excluded_FaceXZ_inp =  np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin"]], dtype=object) #the array of the Rubbish faces coloured as green because they are deemed correct
Cleaned_FaceXZ, Excluded_FaceXZ, indices_to_delete = Delete_clashing_Faces_in_XZandYZ(Cleaned_FaceXZ_inp, Excluded_FaceXZ_inp, FacesXZ)

Cleaned_FaceYZ_inp= np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: y", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin"]], dtype=object)
Excluded_FaceYZ_inp =  np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin"]], dtype=object)
Cleaned_FaceYZ, Excluded_FaceYZ, indices_to_delete = Delete_clashing_Faces_in_XZandYZ(Cleaned_FaceYZ_inp, Excluded_FaceYZ_inp, FacesYZ)

Cleaned_FaceXY_inp= np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Ymax", "Ymin"]], dtype=object)
Excluded_FaceXY_inp =  np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Face Normal", "Face Vertices","Xmax", "Xmin", "Ymax", "Ymin"]], dtype=object)
Cleaned_FaceXY, Excluded_FaceXY, indices_to_delete = Delete_clashing_Faces_in_XYplane(Cleaned_FaceXY_inp, Excluded_FaceXY_inp, FacesXY)
# Using comparison
def get_sign(value):
    if value > 0:
        return 1
    elif value < 0:
        return -1
    else:
        return 0
# Here we match segments and faces. Be carefull two or more faces can be matched two one segment, that's why we need to further process the array: Possible_Correct_Faces
def Assign_Face_to_Segment(Correct_Matching_Segment, Missing_Faces, Cleaned_Faces, Segements):
    #Initialise the arrays, in case in which the for loop does not produce any values for the Possible_Seg_in_Same_AxisNew array
    Possible_Seg_in_Same_AxisNew = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Xmax", "Xmin", "Zmax", "Zmin", "Colour"]], dtype=object)
    for face in Cleaned_Faces[1:]:
     #if face[3] == "05hIZEszH5GBM5fSI338va":
        Possible_Seg_in_Same_Axis = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Xmax", "Xmin", "Zmax", "Zmin", "Colour"]], dtype=object)
        for segment in Segements[1:]:
            AverageLocation_difference_1 = abs(abs(face[7]) - abs(segment[2])) # this is the primarily location difference between a face and a segment, namely in the direction of the face's normal
            Condition1_Axis2 = abs(abs(face[10]) - abs(segment[5]))
            Condition2_Axis3 = abs(abs(face[11]) - abs(segment[6]))
            condition3_Axis2 = 2*segment[3] - face[8]
            condition3_Axis3 = 2*segment[4] - face[9]
            Face_normal = face[12]
            Normal_Y_direction = Face_normal[1]
            if face[2] == "IfcColumn": #If IfcColumn then we modify the condition 3 to say that include the segments that have for example dysegment <= 3*dyface/2. This condition based on the hypothesis that only Column segments can satisfy this condition.
                condition3_Axis2 = 3*face[8]/2 - segment[3]
                if Condition1_Axis2 <= face[8]/2 and Condition2_Axis3 <= face[9]/2 and condition3_Axis2 >=0 and condition3_Axis3 >= 0 : # the "and" in this function denotes that each condition creates a subset of segments, by set that we want to satisfy every condition we reduce the subset of segments to one segment.
                    New_Corr_face_object = np.array([face[0], face[1],face[2], face[3], face[4], face[5], face[6], face[7], face[8], face[9], face[10], face[11], face[12], face[13], face[14], face[15], face[16], face[17], segment[0], segment[1], segment[2], segment[3], segment[4], segment[5], segment[6], segment[7], segment[8], segment[9], segment[10], segment[11]], dtype=object)
                    Possible_Seg_in_Same_Axis = np.vstack([Possible_Seg_in_Same_Axis, New_Corr_face_object]) 
            else:
                if Condition1_Axis2 <= face[8]/2 and Condition2_Axis3 <= face[9]/2 and condition3_Axis2 >=0 and condition3_Axis3 >= 0 : # the "and" in this function denotes that each condition creates a subset of segments, by set that we want to satisfy every condition we reduce the subset of segments to one segment.
                    New_Corr_face_object = np.array([face[0], face[1],face[2], face[3], face[4], face[5], face[6], face[7], face[8], face[9], face[10], face[11], face[12], face[13], face[14], face[15], face[16], face[17], segment[0], segment[1], segment[2], segment[3], segment[4], segment[5], segment[6], segment[7], segment[8], segment[9], segment[10], segment[11]], dtype=object)
                    Possible_Seg_in_Same_Axis = np.vstack([Possible_Seg_in_Same_Axis, New_Corr_face_object]) 
        #Condition 4_0: this is The most crucial condition of all. It sets the distance in Axis 2 in which a correct segment should be located
        #So Possible_Seg_in_Same_Axis array might contain many segments in the direction of Axis 2
        #Here we will limit the segments to maximum 2. This is inevtable becasue the segments located in the same Building element and in the same plane can have differenc less than 0.3 m in axis 2
        #Becasue in our case the maximum wall thickness is 0.3m, to include both the segments of the same Building Element we increse the tolerance of Condition 4_0 to 0.34. But our minimum distnace from wall to wall in the IFC model is 0.595 so if we increase the tolerance to 0.34 a segment might be assigned to two walls. So Condition New: Wall to wall distance/2 < min(Wall_Thickness), but let;s continue to find the results and we will come back to this later
        Possible_Seg_in_Same_AxisNew = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Xmax", "Xmin", "Zmax", "Zmin", "Colour"]], dtype=object)
        i = 1 
        for row in Possible_Seg_in_Same_Axis[1:]:
            Axis1_diff = abs(abs(row[7]) - abs(row[20]))
            if Axis1_diff <= 0.34: #then we are 100% sure that there are only two possible segments for the face under investigation, This value can change according to the nature of the example we investigate.
                New_Corr_row_object = np.array([row[0], row[1],row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20], row[21], row[22], row[23], row[24], row[25], row[26], row[27], row[28], row[29]], dtype=object)
                Possible_Seg_in_Same_AxisNew = np.vstack([Possible_Seg_in_Same_AxisNew, New_Corr_row_object])
                print(Possible_Seg_in_Same_AxisNew)
            else:
                i = i + 1
            #if i == 3:
            #print("the two possible segments do not satsify the condition 4, so the face is labeled as missing")
            #face_colour = [1,0,0] #(Red, Green, Blue)
            #New_Missing_face_object = np.array([row[0], row[1],row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], face_colour ], dtype=object)
            #Missing_Faces = np.vstack([Missing_Faces, New_Missing_face_object])
                    
        #So far we have created a new array Possible Correct segments for one face in the same axis. 
        # Now if the current face has only one row in this array means that the array includes only the Headers and the algorithm did not find any possible segments, so
        if len(Possible_Seg_in_Same_AxisNew) == 1:
            print("Array contains only header row and the algorithm haven't identified any matching segment")
            face_colour = [1,0,0] #(Red, Green, Blue)
            New_Missing_face_object = np.array([face[0], face[1],face[2], face[3], face[4], face[5], face[6], face[7], face[8], face[9], face[10], face[11], face[12], face[13], face[14], face[15], face[16], face[17], face_colour ], dtype=object)
            Missing_Faces = np.vstack([Missing_Faces, New_Missing_face_object])
        #else there is at least one matching segment for the current face 
        # if the length of the Possible_Seg_in_Same_Axis array is 2 means that there is only one correct segment for this face.
        elif len(Possible_Seg_in_Same_AxisNew) == 2:
            Correct_Matching_Segment = np.vstack([Correct_Matching_Segment, Possible_Seg_in_Same_AxisNew[1]])
            
        else:
            twin_faces_array = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin"]], dtype=object)
            #Find twin face of each face
            value_to_compare_IfcType = Possible_Seg_in_Same_AxisNew[1, 2]
            value_to_compare_GUID = Possible_Seg_in_Same_AxisNew[1, 3]
            for face1 in Cleaned_Faces[1:]: #Find the twin face of this face in Cleaned_Faces array
                #In cleaned faces we will always have two faces for each element, as every time we insert only the Clenaed faces for XZ, YZ, and XY planes, so the maximum number of faces that this array can contain for each element is 2.
                #So we add those faces in a new array. If there are two faces then the face under investigation has twin.
                #In this for loop we care only about the faces geometry not segments geometry
                if value_to_compare_IfcType == face1[2] and value_to_compare_GUID == face1[3]:
                    New_twin_face_object = np.array([face1[0], face1[1],face1[2], face1[3], face1[4], face1[5], face1[6], face1[7], face1[8], face1[9], face1[10], face1[11], face1[12], face1[13], face1[14], face1[15], face1[16], face1[17]], dtype=object)
                    twin_faces_array = np.vstack([twin_faces_array, New_twin_face_object])
            #Now we have the twin faces in an array
            #The array Possible_Seg_in_Same_Axis, contains the geometry of the possible segments that we care about, however, this array might contains also segments from the elements of the opposite side of the room, namely segments that we do not care about
            # To avoid considering those different element segments as correct, we hupothesise that the smallest dimension of a room in as-built model is 0.6m so the maximum relative deviation that the two opposite elements can have is 0.3m, 
            #so We create a new Possible_Seg_in_Same_Axis array that should contain only the possible segments for our face.
            
                    
            #Now we are 100%sure that the Possible_Seg_in_Same_AxisNew contain only two possible segemnts for the face under investigation
            #put in a same array the face under investigation and its twin.
            #This procedure is done to be able to model the 6 different cases of clashing between a face and a segment.
            twin_faces_arrayNew = np.empty_like(twin_faces_array)
            twin_faces_arrayNew[0] = twin_faces_array[0]          
            for row1 in twin_faces_array[1:]:
                if face[1] != row1[1]:
                    twin_faces_arrayNew[2] = row1 
                else:
                    twin_faces_arrayNew[1] = row1
            #in this way we will always have in the first row the face that is under investigation
            print(twin_faces_arrayNew)
            #Now we can start Investigating the different cases of 
            if len(Possible_Seg_in_Same_AxisNew) > 1:
                Aver_Locat_seg1 = Possible_Seg_in_Same_AxisNew[1, 20]
                Aver_Locat_seg2 = Possible_Seg_in_Same_AxisNew[2, 20]
                Max_segment = max(Aver_Locat_seg1, Aver_Locat_seg2)
                if Aver_Locat_seg2 == Max_segment:
                    Aver_Locat_seg2 = Aver_Locat_seg1
                    Aver_Locat_seg1 = Max_segment
                #With this way the variable Aver_Locat_seg1 will always represent the Segment 1, namely the segment with the largest average location in Axis 1 between seg 1 and 2.
                #Accordingly, the  Aver_Locat_seg2 variable will always represent the Segment 2, namely the segment with the lowest  average location in Axis 1 between seg 1 and 2.
                    
                Aver_Locat_face1 = twin_faces_arrayNew[1, 7] #This will always be the current face under investigation
                Aver_Locat_face2 = twin_faces_arrayNew[2, 7] #This will always be the twin face of the current face under investigation
                def Assign_segment_to_face(Aver_Locat_seg1, Aver_Locat_seg2, Aver_Locat_face1, Aver_Locat_face2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment ):
                    if Aver_Locat_face1 >= Aver_Locat_face2:
                        #Then the segment that fits to the face under investigation (face XZ 1) is the segment that has the highest value between the segment 1 and segment 2.
                        Matching_Segment = max(Aver_Locat_seg1, Aver_Locat_seg2)
                        for row in Possible_Seg_in_Same_AxisNew[1:]:
                            if Matching_Segment == row[20]:
                                Correct_Matching_Segment = np.vstack([Correct_Matching_Segment, row])
                    else:
                        Matching_Segment = min(Aver_Locat_seg1, Aver_Locat_seg2)
                        for row in Possible_Seg_in_Same_AxisNew[1:]:
                            if Matching_Segment == row[20]:
                                Correct_Matching_Segment = np.vstack([Correct_Matching_Segment, row])
                    return Correct_Matching_Segment
                def Max_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment ):
                    #Then the segment that fits to the face under investigation (face XZ 1) is the segment that has the highest value between the segment 1 and segment 2.
                    Matching_Segment = max(Aver_Locat_seg1, Aver_Locat_seg2)
                    for row in Possible_Seg_in_Same_AxisNew[1:]:
                        if Matching_Segment == row[20]:
                            Correct_Matching_Segment = np.vstack([Correct_Matching_Segment, row])
                    return Correct_Matching_Segment
                
                def Min_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment ):
                    Matching_Segment = min(Aver_Locat_seg1, Aver_Locat_seg2)
                    for row in Possible_Seg_in_Same_AxisNew[1:]:
                        if Matching_Segment == row[20]:
                            Correct_Matching_Segment = np.vstack([Correct_Matching_Segment, row])
                    return Correct_Matching_Segment
                #Case 1: If both seg1, and seg2 average location in Axis 1 is bigger than both face 1 and 2 average locations then:
                if Aver_Locat_seg1 > Aver_Locat_face1 and Aver_Locat_seg2 > Aver_Locat_face1 and Aver_Locat_seg1 > Aver_Locat_face2 and Aver_Locat_seg2 > Aver_Locat_face2:
                    Correct_Matching_Segment = Assign_segment_to_face(Aver_Locat_seg1, Aver_Locat_seg2, Aver_Locat_face1, Aver_Locat_face2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
                #Case 2
                elif Aver_Locat_seg1 < Aver_Locat_face1 and Aver_Locat_seg2 < Aver_Locat_face1 and Aver_Locat_seg1 < Aver_Locat_face2 and Aver_Locat_seg2 < Aver_Locat_face2:
                    Correct_Matching_Segment = Assign_segment_to_face(Aver_Locat_seg1, Aver_Locat_seg2, Aver_Locat_face1, Aver_Locat_face2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
                #Case 3
                elif Aver_Locat_seg1 > Aver_Locat_face1 and Aver_Locat_seg2 < Aver_Locat_face1 and Aver_Locat_seg1 > Aver_Locat_face2 and Aver_Locat_seg2 > Aver_Locat_face2:
                    Correct_Matching_Segment = Max_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
                #Case 4
                elif Aver_Locat_seg1 > Aver_Locat_face1 and Aver_Locat_seg2 > Aver_Locat_face1 and Aver_Locat_seg1 > Aver_Locat_face2 and Aver_Locat_seg2 < Aver_Locat_face2:
                    Correct_Matching_Segment = Min_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
                #Case 5
                elif Aver_Locat_seg1 < Aver_Locat_face1 and Aver_Locat_seg2 < Aver_Locat_face1 and Aver_Locat_seg1 > Aver_Locat_face2 and Aver_Locat_seg2 < Aver_Locat_face2:
                    Correct_Matching_Segment = Max_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
                #Case 6
                elif Aver_Locat_seg1 > Aver_Locat_face1 and Aver_Locat_seg2 < Aver_Locat_face1 and Aver_Locat_seg1 < Aver_Locat_face2 and Aver_Locat_seg2 < Aver_Locat_face2:
                    Correct_Matching_Segment = Min_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
                #Case 7 
                elif Aver_Locat_seg1 < Aver_Locat_face1 and Aver_Locat_seg2 < Aver_Locat_face1 and Aver_Locat_seg1 > Aver_Locat_face2 and Aver_Locat_seg2 > Aver_Locat_face2:
                    Correct_Matching_Segment = Max_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
                #Case 8
                elif Aver_Locat_seg1 < Aver_Locat_face2 and Aver_Locat_seg2 < Aver_Locat_face2 and Aver_Locat_seg1 > Aver_Locat_face1 and Aver_Locat_seg2 > Aver_Locat_face1:
                    Correct_Matching_Segment = Min_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
                #Case 9
                elif Aver_Locat_seg1 > Aver_Locat_face2 and Aver_Locat_seg2 < Aver_Locat_face2 and Aver_Locat_seg1 > Aver_Locat_face1 and Aver_Locat_seg2 < Aver_Locat_face1 and Aver_Locat_face1 > Aver_Locat_face2:
                    Correct_Matching_Segment = Max_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
                 #Case 10
                elif Aver_Locat_seg1 > Aver_Locat_face2 and Aver_Locat_seg2 < Aver_Locat_face2 and Aver_Locat_seg1 > Aver_Locat_face1 and Aver_Locat_seg2 < Aver_Locat_face1 and Aver_Locat_face1 < Aver_Locat_face2:
                    Correct_Matching_Segment = Min_Segment(Aver_Locat_seg1, Aver_Locat_seg2, Possible_Seg_in_Same_AxisNew, Correct_Matching_Segment)
    return Correct_Matching_Segment, Missing_Faces, Possible_Seg_in_Same_AxisNew

# Find the Possible Correct faces that might be match with each Face
Missing_FacesXZ = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin", "Colour" ]], dtype=object)
Correct_Matching_SegmentsXZ = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Xmax", "Xmin", "Zmax", "Zmin", "Colour" ]], dtype=object)
Correct_Matching_SegmentsXZ, Missing_FacesXZ, Possible_Seg_in_Same_AxisXZ = Assign_Face_to_Segment(Correct_Matching_SegmentsXZ, Missing_FacesXZ, Cleaned_FaceXZ, Segement_XZ)
print(Correct_Matching_SegmentsXZ)

Missing_FacesYZ = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin", "Colour" ]], dtype=object)
Correct_Matching_SegmentsYZ = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Ymax", "Ymin", "Zmax", "Zmin", "Colour" ]], dtype=object)
Correct_Matching_SegmentsYZ, Missing_FacesYZ, Possible_Seg_in_Same_AxisYZ = Assign_Face_to_Segment(Correct_Matching_SegmentsYZ, Missing_FacesYZ, Cleaned_FaceYZ, Segement_YZ)
print(Correct_Matching_SegmentsYZ)

Missing_FacesXY = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Ymax", "Ymin", "Colour" ]], dtype=object)
Correct_Matching_SegmentsXY = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Ymax", "Ymin", "Segment ID", "Segment Name", "Segment Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Xmax", "Xmin", "Ymax", "Ymin", "Colour" ]], dtype=object)
Correct_Matching_SegmentsXY, Missing_FacesXY, Possible_Seg_in_Same_AxisXY  = Assign_Face_to_Segment(Correct_Matching_SegmentsXY, Missing_FacesXY, Cleaned_FaceXY, Segement_XY)
print(Correct_Matching_SegmentsXY)
# We need to take into account also for the Rubbish faces or for the face that have not any segment any near these faces will be absent from the as-built model
#Compare a face with a segment
# Finally we can also load the corresponding segment of each face in the mesh lab for better visualisation
def Colour_allocation(Deviation_in_axis):
    if abs(Deviation_in_axis) <= 0.05:
        face_colour = [0,1,0]# Green colour (RGB)=(Red, Green, Blue)
    else:
        face_colour = [1,0,0]
    return face_colour
def Colour_allocation_OE(Deviation_in_axis):
    if abs(Deviation_in_axis) <= 0.05:
        face_colour = [0,153,0]# Green colour (RGB)=(Red, Green, Blue)
    else:
        face_colour = [153,0,0]
    return face_colour
#for plane XZ, see the columns with indices 7,8,9 of the Correct_Matching_SegmentsXZ array. 
#The corresponding Axis 1: Y-axis, Axis 2: X-axis, Axis 3: Z-axis, represented in Index 7,8,9 respectively. the same for the rest of planes YZ, and XY
#for every plane Axis 1, always represent the coordinate which all point of the face or segment have in common, in case of XZ plane, y-coordinate
#For every plane Axis 1, and Axis 2,  represent the dimensions of the Plane. In case of XZ plane: Axis 2: Dx, Axis 3: Dz
#In the below function, we calculate the deviations in the X, Y, Z axes between a face and its corresponding segment. 
def Deviations_in_X_Y_Z_axes(Correct_Matching_Segments, Final_Face_Segment):
    New_Missing_face_object = np.array([None, None, None, None, None], dtype=object) #in this array all the face colours for each of the five comparisons where stored. In the end we will check if at least one of them is red, then  the entire face is coloured red, becasue one of the five deviations is exceeding tolerance
    for row in Correct_Matching_Segments[1:]:
     #if row[3] == "05hIZEszH5GBM5fSI338va":
        #for axis 1
        if row[20] >= row[7]:
            Deviation_in_1_axis = abs(row[20] - row[7]) #Deviation in positive corresponding axis, the sign of the deviation denoted the direction of the deviation
        else:
            Deviation_in_1_axis = -abs(row[20] - row[7]) #Deviation in negative corresponding axis, the sign of the deviation denoted the direction of the deviation
        face_colour1 = Colour_allocation(Deviation_in_1_axis)
        New_Missing_face_object[0] = face_colour1 
        #for axis 2 Now we compare with Max and Min coordinates of the axis
        #For Axis 2- MAX
        #The logic is that we take the Max coordinate in the  axis 2 of the face, and we compare with the max coordinate in the axis 2 of the segment, 
        #then we  compare the face with the segment, but be carefull, we consider the face as static, and the segment as moving, so 
        # we say if Segment max is bigger than segment max, then the face has been moved in the positive direction of Axis 2 by Deviation_in_Max_in_axis2.
        if row[25] >= row[14]:  #Example faceXZ: XfaceMax=row[14]>= XSegMax=row[25], 
            Deviation_in_Max_in_axis2 = abs(row[25] - row[14]) - 0.015 #Deviation in positive corresponding axis, we account also for The segmentation algorithm errors
        else:
            Deviation_in_Max_in_axis2 = -abs(row[25] - row[14]) + 0.015 #deviation in negative corresponding axis, we account also for The segmentation algorithm errors
        face_colour2 = Colour_allocation(Deviation_in_Max_in_axis2)
        New_Missing_face_object[1] = face_colour2
        #For Axis 2- MIN
        if row[26] >= row[15]:  #Example faceXZ: XfaceMin=row[15]>= XSegMin=row[26]
            Deviation_in_Min_in_axis2 = abs(row[26] - row[15]) - 0.015 #Deviation in positive corresponding axis, we account also for The segmentation algorithm errors
        else:
            Deviation_in_Min_in_axis2 = -abs(row[26] - row[15])  + 0.015 #deviation in negative corresponding axis, we account also for The segmentation algorithm errors
        face_colour3 = Colour_allocation(Deviation_in_Min_in_axis2)
        New_Missing_face_object[2] = face_colour3
        
        #For Axis 3- MAX
        if row[27] >= row[16]:
            Deviation_in_Max_in_axis3 = abs(row[27] - row[16]) #Deviation in positive corresponding axis
        else:
            Deviation_in_Max_in_axis3 = -abs(row[27] - row[16]) #deviation in negative corresponding axis
        face_colour4 = Colour_allocation(Deviation_in_Max_in_axis3)
        New_Missing_face_object[3] = face_colour4
        #For Axis 3- MIN
        if row[28] >= row[17]:
            Deviation_in_Min_in_axis3 = abs(row[28] - row[17]) #Deviation in positive corresponding axis
        else:
            Deviation_in_Min_in_axis3 = -abs(row[28] - row[17]) #deviation in negative corresponding axis
        face_colour5 = Colour_allocation(Deviation_in_Min_in_axis3)
        New_Missing_face_object[4] = face_colour5
        face_colour_final = [0,1,0] # We set the final face colour as green by default. Then we are looking if there is any red colour produced from the five types of deviation
        for i in range(len(New_Missing_face_object)):
            #if at least one of the colours are red then at least one within the five types of deviations is exceed the tolerance
            if New_Missing_face_object[i] == [1,0,0]:
                face_colour_final = New_Missing_face_object[i] #then the face is coloured Red.
                
        #when the for finish, 
        New_Final_face_object = np.array([row[0], row[1],row[2], row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20], row[21], row[22], row[23], row[24], row[25], row[26], row[27], row[28], face_colour_final, Deviation_in_1_axis, Deviation_in_Max_in_axis2, Deviation_in_Min_in_axis2, Deviation_in_Max_in_axis3, Deviation_in_Min_in_axis3 ], dtype=object)
        Final_Face_Segment = np.vstack([Final_Face_Segment, New_Final_face_object])
    return Final_Face_Segment
#Final_Face_Segment arrays are actually the Correct_Matching_Segment arrays but with info about Comparison results
Final_Face_SegmentXZ = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Xmax", "Xmin", "Zmax", "Zmin", "Colour", "Deviation in Y-axis", "Deviation of face Xmax", "Deviation of face Xmin",  "Deviation of face Zmax", "Deviation of face Zmin" ]], dtype=object) # the axes in these arrays refer to the GCS
Final_Face_SegmentXZ = Deviations_in_X_Y_Z_axes(Correct_Matching_SegmentsXZ, Final_Face_SegmentXZ)
        
Final_Face_SegmentYZ = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Ymax", "Ymin", "Zmax", "Zmin", "Colour", "Deviation in X-axis", "Deviation of face Ymax", "Deviation of face Ymin",  "Deviation of face Zmax", "Deviation of face Zmin" ]], dtype=object) #35 Columns, 34 indices
Final_Face_SegmentYZ = Deviations_in_X_Y_Z_axes(Correct_Matching_SegmentsYZ, Final_Face_SegmentYZ)    

Final_Face_SegmentXY = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Ymax", "Ymin", "Segment ID", "Segment Name", "Segment Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Xmax", "Xmin", "Ymax", "Ymin", "Colour", "Deviation in Z-axis", "Deviation of face Xmax", "Deviation of face Xmin",  "Deviation of face Ymax", "Deviation of face Ymin" ]], dtype=object)
Final_Face_SegmentXY = Deviations_in_X_Y_Z_axes(Correct_Matching_SegmentsXY, Final_Face_SegmentXY) 
print(Final_Face_SegmentXZ)

#COMPARISON OF IFCOPENINGELEMENTS
def Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_minus, Axis2_max, Axis3_minus, Axis3_max):
    if face_OE[0]=="Face XZ":
        #For OpenElement in XZ axis: Axis 1 = Y-axis, Axis 2 = X-axis, and Axis 3 = Z-axis
        selected_indices = np.where(
            (points[:, 0] >= Axis2_minus) & (points[:, 0] <= Axis2_max) &
            (points[:, 1] >= Axis1_minus) & (points[:, 1] <= Axis1_max) &
            (points[:, 2] >= Axis3_minus) & (points[:, 2] <= Axis3_max))[0]
        selected_points = points[selected_indices]
    elif face_OE[0]=="Face YZ":
        #For OpenElement in YZ axis: Axis 1 = X-axis, Axis 2 = Y-axis, and Axis 3 = Z-axis
        selected_indices = np.where(
            (points[:, 0] >= Axis1_minus) & (points[:, 0] <= Axis1_max) &
            (points[:, 1] >= Axis2_minus) & (points[:, 1] <= Axis2_max) &
            (points[:, 2] >= Axis3_minus) & (points[:, 2] <= Axis3_max))[0]
        selected_points = points[selected_indices]
    return selected_points
def Segemnt_Coordinates_in_BB(points_Bound_Box, face_OE):  
    # we Can have two cases: 1) face on XZ plane and 2) face on YZ plane
    #First we place in seperate arrays the coordinates x, y, z of the selected point are contained in the bounding box 
    if points_Bound_Box.size == 0:
       Axis1_seg_max = 0
       Axis2_seg_max = 0
       Axis3_seg_max = 0
       Axis1_seg_min = 0
       Axis2_seg_min = 0
       Axis3_seg_min = 0
    else:
        Seg_x = points_Bound_Box[:, 0] #This is slicing operation selects all rows (:) and only the elements in the first column (0)
        Seg_y = points_Bound_Box[:, 1]
        Seg_z = points_Bound_Box[:, 2]

        if face_OE[0]=="Face XZ": #If the face is located in XZ plane the Axis 1, 2, 3 is the Y, X, Z axes respectively
            #Max values of segment points in BB1 
            Axis1_seg_max = max(Seg_y)
            Axis2_seg_max = max(Seg_x)
            Axis3_seg_max = max(Seg_z)
            #Min values of segment points in BB1 
            Axis1_seg_min = min(Seg_y)
            Axis2_seg_min = min(Seg_x)
            Axis3_seg_min = min(Seg_z)
        elif face_OE[0]=="Face YZ": #If the face is located in YZ plane the Axis 1, 2, 3 is the X, Y, Z axes respectively
            #Max values of segment points in BB1 
            Axis1_seg_max = max(Seg_x)
            Axis2_seg_max = max(Seg_y)
            Axis3_seg_max = max(Seg_z)
            #Max values of segment points in BB1 
            Axis1_seg_min = min(Seg_x)
            Axis2_seg_min = min(Seg_y)
            Axis3_seg_min = min(Seg_z)
    return Axis1_seg_max, Axis2_seg_max, Axis3_seg_max, Axis1_seg_min, Axis2_seg_min, Axis3_seg_min
#In correct matching Segment array there are only the wall faces that are already existing, therefore, if a Opening element Face does not have any match, it is automatically classified to missing elements
def OpeningElement_Comparison(Correct_Matching_Segment, Corrrect_Mat_Seg_2nd_axis, Faces_OE,  Missing_Faces, Final_Face_Segment, input_Segment_directory_2):
    # Filter rows where "IfcType" is "IfcWall"
    filtered_rows = Correct_Matching_Segment[1:][Correct_Matching_Segment[1:, 2] == "IfcWall"]
    # Combine the header with the filtered rows
    IfcWall_rows = np.vstack((Correct_Matching_Segment[0], filtered_rows))
    for face_OE in Faces_OE[1:]: 
     #if face_OE[3] == "1mNm2D_bLEcw7ElR7oHOnv":
        i = 1 #We initialise this variable to avoid the algorithm goes into the if again after has found the wall tha matches with the opening Element
        for face_Wall in IfcWall_rows[1:]:
            if re.search(patternIfcWall, str(face_Wall[2])) and face_OE[0] == face_Wall[0] and face_OE[1] == face_Wall[1] and face_Wall[5] < face_OE[5] < face_Wall[5] + 1: #Here we find the Wall face that is assigned to the current OpenElement face by checking if the OpeningElement's Local ID is based on Wall's face Local ID, Ex: Wallface ID = 1, then all its opening Elements should have IDs, 1.1, 1.2, 1.3, etc
                New_Colour_face_object = np.array([None, None, None, None, None], dtype=object) #in this array all the face colours for each of the five comparisons where stored. In the end we will check if at least one of them is red, then  the entire face is coloured red, becasue one of the five deviations is exceeding tolerance
                New_FaceOE = None
                face_colour = [0,153,0]
                i = i + 1 # by increasing this variable to 2 we stop the algorithm to insert in this if command for second time because it will write down the face as missign face.
                New_FaceOE = np.array([face_OE[0], face_OE[1],face_OE[2], face_OE[3], face_OE[4], face_OE[5], face_OE[6], face_OE[7], face_OE[8], face_OE[9], face_OE[10], face_OE[11], face_OE[12], face_OE[13], face_OE[14], face_OE[15], face_OE[16], face_OE[17], face_Wall[18], face_Wall[19], face_Wall[20], face_Wall[21], face_Wall[22], face_Wall[23], face_Wall[24], face_Wall[25], face_Wall[26], face_Wall[27], face_Wall[28], face_Wall[29]], dtype=object)
                #Now we have in the array: New_FaceOE_Object the info about the face OE under investigation and the info about the segement that is assigned to it
                #But we need the point cloud file of the current Segment, for that we use the New ply files that we created in directory: output_Segment_directory, which are named with Unique way
                #Construct the filename for the desirable segment PLY file
                segment_filename = f"{New_FaceOE[19]}_{New_FaceOE[18]}.ply"
                for filename in os.listdir(input_Segment_directory_2):    
                    if filename == segment_filename:
                        # Load .ply file
                        segment_filepath = os.path.join(input_Segment_directory_2, filename)
                        # Load the corresponding point cloud file
                        segment_ply = o3d.io.read_point_cloud(segment_filepath)
                        # Get points as numpy array
                        points = np.asarray(segment_ply.points)
                #extract the four vertices of the OE face
                face_vertices_array = face_OE[13]
                # Extract x, y, and z coordinates
                x = face_vertices_array[:, 0] #This is slicing operation selects all rows (:) and only the elements in the first column (0)
                y = face_vertices_array[:, 1]
                z = face_vertices_array[:, 2]
                # Calculate differences
                dx = np.max(x) - np.min(x)
                dy = np.max(y) - np.min(y)
                dz = np.max(z) - np.min(z)
                #Create a bounding box within the 8 vertices of the IfcOpeningElement
                # Find min and max values of the Bounding Box for Axis 1
                Axis1_minus = face_OE[7] - 0.1
                Axis1_max = face_OE[7] + 0.1
                # Find min and max values of the Bounding Box for Axis 2
                Axis2_minus = face_OE[15]
                Axis2_max = face_OE[14]
                # Find min and max values of the Bounding Box for Axis 3
                Axis3_max = face_OE[16]
                Axis3_minus = face_OE[17]
                # Select points within the cube boundaries
                selected_points = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_minus, Axis2_max, Axis3_minus, Axis3_max)
                # The Axis 1,2,3 can change as the Opening Element can located in XZ and YZ
                #Now we have created the bounding boxes, and we have selected the points that are contained inside.
                #so we are ready for the Comparison

                #we create Two new bounding Boxes
                #From now on: Edge_Max --> Edge2 max, and Edge_Min --> Edge1 min
                #For Deviations in Axis 2, and from the edges and inwards, we should create two more Bounding boxes as subsets of the first bounding box
                Axis2_Average = Axis2_minus + (Axis2_max - Axis2_minus)/2
                Axis2_Edge1_max = Axis2_Average - 0.1
                Axis2_Edge1_min = Axis2_minus
                Axis2_Edge2_max = Axis2_max
                Axis2_Edge2_min = Axis2_Average + 0.1
                #Now we have the boundaries of the two new bounding boxes
                #25/05/2024: we have realised that we need to split the previous two BBs into two BBs each, so now we have 4 BBs in order to account also for Axis 3

                Axis3_Average = Axis3_minus + (Axis3_max - Axis3_minus)/2
                Axis3_Edge_Z1_max = Axis3_max
                Axis3_Edge_Z1_min = Axis3_Average + 0.1
                Axis3_Edge_Z2_max = Axis3_Average - 0.1
                Axis3_Edge_Z2_min = Axis3_minus
                #Now we have the boundaries of the four new bounding boxes in all axes
                points_Bound_Box_1_1 = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge1_min, Axis2_Edge1_max, Axis3_Edge_Z1_min, Axis3_Edge_Z1_max)
                points_Bound_Box_1_2 = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge1_min, Axis2_Edge1_max, Axis3_Edge_Z2_min, Axis3_Edge_Z2_max)
                
                points_Bound_Box_2_1 = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge2_min, Axis2_Edge2_max, Axis3_Edge_Z1_min, Axis3_Edge_Z1_max)
                points_Bound_Box_2_2 = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge2_min, Axis2_Edge2_max, Axis3_Edge_Z2_min, Axis3_Edge_Z2_max)
                #After we have seperated the two boxes, we managed to isolate the points that are related to the EdgeMin, namely face[15],from Bounding Box 1,  and Edge Max, namely face[14], from Bounding Box 2
                
                #For bounding Box 1.1 we calculate two deviations.
                if points_Bound_Box_1_1.size > 0: #this means that array is not empty, then it might have inward or outward deviation 
                    # The Bounding Box 1.1 is extended in Axis 2 between Axis2_Edge1_min and Axis2_Edge1_max,
                    #and in Axis 3 between Edge_Z1_min and Edge_Z1_max
                    Axis1_seg1_max, Axis2_seg1_max, Axis3_seg1_max, Axis1_seg1_min, Axis2_seg1_min, Axis3_seg1_min = Segemnt_Coordinates_in_BB(points_Bound_Box_1_1, face_OE)
                    #Now we have all the values that represent the points of the segment that are included inside each bounding box 1.1. Apparently we are not going to use all values.
                    #We want to maintain the values of these variables so we have to be very careful to not replace their names
                    
                    #IMPORTANT: We set maximum deviation that the algorithm can identified as:
                    Max_deviation = Axis2_Edge1_max - Axis2_Edge1_min
                    #2) FOR AXIS 2
                    #For Deviations in Axis 2, and from the edge and inwards:
                    #We find the lengths, first of Face in axis 3, and of the segment in Axis 3 for the points included in BB 1.1
                    Axis_3_seg_Diff = Axis3_seg1_max - Axis3_seg1_min #From initial points inside the BB 1.1
                    Axis_3_face_Diff = Axis3_Edge_Z1_max - Axis3_Edge_Z1_min
                    Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                    Deviation_in_Min_in_axis2 = 0
                    Axis2_Edge1_min_new = Axis2_Edge1_min
                    Axis2_Edge1_min_initial = Axis2_Edge1_min
                    points_Bound_Box_1_1_new = points_Bound_Box_1_1.copy()
                    if Face_seg_diff <= 0.05: #This means that the points are extended in the entire length of the face in BB 1.1 in Axis 2 so, we are usre that we have deviation in Z axis, as the entire wall has moved as a whole
                        Step_sum = 0
                        while Face_seg_diff <= 0.05:# We are looking for inward deviation
                            #in this loop we add 0.005 m, namely 1 cm every time from the value of the left boundary of the BB 1.1. With this way we lengthen the size of the BB 1.1 until the above condition not satisfied
                            #We are in BB1.1 so we add 0.005 in the variable: Axis2_Edge1_min
                            Axis2_Edge1_min_new = Axis2_Edge1_min_new + 0.01 #Be careful, We are in the case of Axis 3 ,so we substract 0.01 from the Axis3_Edge_Z1_max
                            #Increase the step by 0.005 to make the algorithm stop if there is no points at all
                            Step_sum = Step_sum + 0.01
                            #then we recalculate the Bounding box by incorporating this adjustment
                            points_Bound_Box_1_1_new = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge1_min_new, Axis2_Edge1_max, Axis3_Edge_Z1_min, Axis3_Edge_Z1_max)
                            Axis1_seg1_max_Dev2, Axis2_seg1_max_Dev2, Axis3_seg1_max_Dev2, Axis1_seg1_min_Dev2, Axis2_seg1_min_Dev2, Axis3_seg1_min_Dev2 = Segemnt_Coordinates_in_BB(points_Bound_Box_1_1_new, face_OE)
                            
                            Axis_3_seg_Diff = Axis3_seg1_max_Dev2 - Axis3_seg1_min_Dev2 # We expect that from a point and then, the length of the segment will change
                            #Be careful, Axis_3_face_Diff this variable remains stable during the loop, it doesn't change
                            Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                            if Step_sum >= Max_deviation:
                                Face_seg_diff = 10
                                Axis2_Edge1_min_new = 0

                    else:
                        Step_sum = 0
                        while Face_seg_diff >= 0.05:# We are looking for outward deviation
                            #in this loop we substract 0.005 m, namely 1 cm every time from the left boundary of the BB 1.1. With this way we reduce the size of the BB 1.1 until the above condition not satisfied
                            #We are in BB1.1 so we substract 0.005 in the variable: Axis2_Edge1_min
                            Axis2_Edge1_min_new = Axis2_Edge1_min_new - 0.005 
                            #Increase the step by 0.005 to make the algorithm stop if there is no points at all
                            Step_sum = Step_sum + 0.005
                            #then we recalculate the Bounding box by incorporating this adjustment
                            points_Bound_Box_1_1_new = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge1_min_new, Axis2_Edge1_max, Axis3_Edge_Z1_min, Axis3_Edge_Z1_max)
                            Axis1_seg1_max_Dev2, Axis2_seg1_max_Dev2, Axis3_seg1_max_Dev2, Axis1_seg1_min_Dev2, Axis2_seg1_min_Dev2, Axis3_seg1_min_Dev2 = Segemnt_Coordinates_in_BB(points_Bound_Box_1_1_new, face_OE)
                            
                            Axis_3_seg_Diff = Axis3_seg1_max_Dev2 - Axis3_seg1_min_Dev2 # We expect that from a point and then, the length of the segment will change
                            #Be careful, Axis_3_face_Diff this variable remains stable during the loop, it doesn't change
                            Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                            if Step_sum >= Max_deviation:
                                Face_seg_diff = 0
                                Axis2_Edge1_min_new = 0
                    Deviation_in_Min_in_axis2 = Axis2_Edge1_min_new - Axis2_Edge1_min_initial  # if Deviation>0 is extended in positive corresponding axis, because always Axis2_Edge1_min > Axis2_Edge1_min_initial, always we have Deviation_in_Min_in_axis2 > 0, if Deviation < 0 is the opposite
                    face_colour = Colour_allocation_OE(Deviation_in_Min_in_axis2)
                    New_Colour_face_object[1] = face_colour
                    print("Deviation Axis 2 min checked")
                    #1) FOR AXIS 3 
                    #For Deviations in Axis 3, and from the edges and inwards:
                    #We find the lengths of Face in axis 2, and of the segment in Axis 2 for the points included in BB 1.1
                    Axis_2_seg_Diff = Axis2_seg1_max - Axis2_seg1_min
                    Axis_2_face_Diff = Axis2_Edge1_max - Axis2_Edge1_min
                    Face_seg_diff = abs(Axis_2_seg_Diff - Axis_2_face_Diff)
                    Deviation_in_Max_in_axis3 = 0
                    Axis3_Edge_Z1_max_new = Axis3_Edge_Z1_max #we save the variable Axis3_Edge_Z1_max in a new variable because we want to maintain the value of Axis3_Edge_Z1_max in ordert o use it again in the next part, in Axis 2 part.
                    Axis3_Edge_Z1_max_initial = Axis3_Edge_Z1_max
                    points_Bound_Box_1_1_new = points_Bound_Box_1_1.copy()
                    if Face_seg_diff <= 0.05: #This means that the points are extended in the entire length of the face in BB 1.1 in Axis 2 so, we are usre that we have deviation in Z axis, as the entire wall has moved as a whole
                        Step_sum = 0
                        while Face_seg_diff <= 0.05:# We are looking for inward deviation
                            #in this loop we substract 0.005 m, namely 1 cm every time from the upper boundary of our BB 1.1. With this way we deduct the size of the BB until the above condition not satisfied
                            #We are in BB1 so we substract 0.005 in the variable: Axis3_Edge_Z1_max
                            Axis3_Edge_Z1_max_new = Axis3_Edge_Z1_max_new - 0.005 #Be careful, We are in the case of Axis 3 ,so we substract 0.01 from the Axis3_Edge_Z1_max_new
                            #Increase the step by 0.005 to make the algorithm stop if there is no points at all
                            Step_sum = Step_sum + 0.005
                            #then we recalculate the Bounding box by incorporating this adjustment
                            points_Bound_Box_1_1_new = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge1_min, Axis2_Edge1_max, Axis3_Edge_Z1_min, Axis3_Edge_Z1_max_new)
                            Axis1_seg1_max_Dev3, Axis2_seg1_max_Dev3, Axis3_seg1_max_Dev3, Axis1_seg1_min_Dev3, Axis2_seg1_min_Dev3, Axis3_seg1_min_Dev3 = Segemnt_Coordinates_in_BB(points_Bound_Box_1_1_new, face_OE)
                            
                            Axis_2_seg_Diff = Axis2_seg1_max_Dev3 - Axis2_seg1_min_Dev3 # We expect that from a point and then, the length of the segment will change
                            #Be careful, Axis_2_face_Diff: this variable remains stable during the loop, it doesn't change
                            Face_seg_diff = abs(Axis_2_seg_Diff - Axis_2_face_Diff)
                            if Step_sum >= Max_deviation:
                                Face_seg_diff = 10
                                Axis3_Edge_Z1_max_new = 0
                    else:
                        Step_sum = 0
                        while Face_seg_diff >= 0.05:# We are looking for outward deviation
                            #in this loop we add 0.005 m, namely 1 cm every time in the upper boundary of our BB 1.1. With this way we extend the size of the BB until the above condition not satisfied
                            #We are in BB1 so we add 0.005 in the value of the variable: Axis3_Edge_Z1_max
                            Axis3_Edge_Z1_max_new = Axis3_Edge_Z1_max_new + 0.005 #Be careful, We are in the case of Axis 3 ,so we substract 0.01 from the Axis3_Edge_Z1_max_new
                            #Increase the step by 0.005 to make the algorithm stop if there is no points at all
                            Step_sum = Step_sum + 0.005
                            #then we recalculate the Bounding box by incorporating this adjustment
                            points_Bound_Box_1_1_new = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge1_min, Axis2_Edge1_max, Axis3_Edge_Z1_min, Axis3_Edge_Z1_max_new)
                            Axis1_seg1_max_Dev3, Axis2_seg1_max_Dev3, Axis3_seg1_max_Dev3, Axis1_seg1_min_Dev3, Axis2_seg1_min_Dev3, Axis3_seg1_min_Dev3 = Segemnt_Coordinates_in_BB(points_Bound_Box_1_1_new, face_OE)
                            
                            Axis_2_seg_Diff = Axis2_seg1_max_Dev3 - Axis2_seg1_min_Dev3 # We expect that from a point and then, the length of the segment will change
                            #Be careful, Axis_2_face_Diff: this variable remains stable during the loop, it doesn't change
                            Face_seg_diff = abs(Axis_2_seg_Diff - Axis_2_face_Diff)
                            #Make the Algorithm stop some time.
                            if Step_sum >= Max_deviation:
                                Face_seg_diff = 0
                                Axis3_Edge_Z1_max_new = 0
                    Deviation_in_Max_in_axis3 = Axis3_Edge_Z1_max_new - Axis3_Edge_Z1_max_initial  #Deviation is extended in negative corresponding axis, because Axis3_Edge_Z1_max < Axis3_Edge_Z1_max_initial then always Deviation_in_Max_in_axis3 < 0, as we desire
                    face_colour = Colour_allocation_OE(Deviation_in_Max_in_axis3)
                    New_Colour_face_object[2] = face_colour
                    print("Deviation Axis 3 max checked")
                else: #This means that the array is empty so there are no inward Deviations, in this case we check for outward deviations, if there is none then the OE edge is considered correct
                    Axis1_seg1_max, Axis2_seg1_max, Axis3_seg1_max, Axis1_seg1_min, Axis2_seg1_min, Axis3_seg1_min = Segemnt_Coordinates_in_BB(points_Bound_Box_1_1, face_OE)
                    #for Axis 2 
                    #Again we need to check for outward deviation in Axis 2
                    #For Deviations in Axis 2, and from the edge and inwards:
                    #We find the lengths, first of Face in axis 3, and of the segment in Axis 3 for the points included in BB 1.1
                    Axis_3_seg_Diff = Axis3_seg1_max - Axis3_seg1_min #From initial points inside the BB 1.1
                    Axis_3_face_Diff = Axis3_Edge_Z1_max - Axis3_Edge_Z1_min
                    Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                    Deviation_in_Min_in_axis2 = 0
                    Axis2_Edge1_min_new = Axis2_Edge1_min
                    Axis2_Edge1_min_initial = Axis2_Edge1_min
                    points_Bound_Box_1_1_new = points_Bound_Box_1_1.copy()
                    if Face_seg_diff >= 0.05: #This means that the points are extended in the entire length of the face in BB 1.1 in Axis 2 so, we are sure that we have deviation in Z axis, as the entire wall has moved as a whole
                        Step_sum = 0
                        while Face_seg_diff >= 0.05:# We are looking for outward deviation
                            #in this loop we substract 0.005 m, namely 1 cm every time from the left boundary of the BB 1.1. With this way we reduce the size of the BB 1.1 until the above condition not satisfied
                            #We are in BB1.1 so we substract 0.005 in the variable: Axis2_Edge1_min
                            Axis2_Edge1_min_new = Axis2_Edge1_min_new - 0.005
                            #Increase the step by 0.005 to make the algorithm stop if there is no points at all
                            Step_sum = Step_sum + 0.005
                            #then we recalculate the Bounding box by incorporating this adjustment
                            points_Bound_Box_1_1_new = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge1_min_new, Axis2_Edge1_max, Axis3_Edge_Z1_min, Axis3_Edge_Z1_max)
                            Axis1_seg1_max_Dev2, Axis2_seg1_max_Dev2, Axis3_seg1_max_Dev2, Axis1_seg1_min_Dev2, Axis2_seg1_min_Dev2, Axis3_seg1_min_Dev2 = Segemnt_Coordinates_in_BB(points_Bound_Box_1_1_new, face_OE)
                            
                            Axis_3_seg_Diff = Axis3_seg1_max_Dev2 - Axis3_seg1_min_Dev2 # We expect that from a point and then, the length of the segment will change
                            #Be careful, Axis_3_face_Diff this variable remains stable during the loop, it doesn't change
                            Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                            #Make the Algorithm stop some time.
                            if Step_sum >= Max_deviation:
                                Face_seg_diff = 0
                                Axis2_Edge1_min_new = 0
                    Deviation_in_Min_in_axis2 = Axis2_Edge1_min_new - Axis2_Edge1_min_initial  # if Deviation>0 is extended in positive corresponding axis, because always Axis2_Edge1_min > Axis2_Edge1_min_initial, always we have Deviation_in_Min_in_axis2 > 0, if Deviation < 0 is the opposite
                    face_colour = Colour_allocation_OE(Deviation_in_Min_in_axis2)
                    New_Colour_face_object[1] = face_colour
                    print("Deviation Axis 2 min checked")
                    #For Axis 3
                    #Again we need to check for outward deviation in Axis 3
                    #For Deviations in Axis 3, and from the edges and inwards:
                    #We find the lengths of Face in axis 2, and of the segment in Axis 2 for the points included in BB 1.1
                    Axis_2_seg_Diff = Axis2_seg1_max - Axis2_seg1_min
                    Axis_2_face_Diff = Axis2_Edge1_max - Axis2_Edge1_min
                    Face_seg_diff = abs(Axis_2_seg_Diff - Axis_2_face_Diff)
                    Deviation_in_Max_in_axis3 = 0
                    Axis3_Edge_Z1_max_new = Axis3_Edge_Z1_max #we save the variable Axis3_Edge_Z1_max in a new variable because we want to maintain the value of Axis3_Edge_Z1_max in ordert o use it again in the next part, in Axis 2 part.
                    Axis3_Edge_Z1_max_initial = Axis3_Edge_Z1_max
                    points_Bound_Box_1_1_new = points_Bound_Box_1_1.copy()
                    if Face_seg_diff >= 0.05: 
                         Step_sum = 0
                         while Face_seg_diff >= 0.05:# We are looking for outward deviation
                            #in this loop we add 0.005 m, namely 1 cm every time in the upper boundary of our BB 1.1. With this way we extend the size of the BB until the above condition not satisfied
                            #We are in BB1 so we add 0.005 in the value of the variable: Axis3_Edge_Z1_max
                            Axis3_Edge_Z1_max_new = Axis3_Edge_Z1_max_new + 0.005 #Be careful, We are in the case of Axis 3 ,so we substract 0.01 from the Axis3_Edge_Z1_max_new
                            #Increase the step by 0.005 to make the algorithm stop if there is no points at all
                            Step_sum = Step_sum + 0.005
                            #then we recalculate the Bounding box by incorporating this adjustment
                            points_Bound_Box_1_1_new = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge1_min, Axis2_Edge1_max, Axis3_Edge_Z1_min, Axis3_Edge_Z1_max_new)
                            Axis1_seg1_max_Dev3, Axis2_seg1_max_Dev3, Axis3_seg1_max_Dev3, Axis1_seg1_min_Dev3, Axis2_seg1_min_Dev3, Axis3_seg1_min_Dev3 = Segemnt_Coordinates_in_BB(points_Bound_Box_1_1_new, face_OE)
                            
                            Axis_2_seg_Diff = Axis2_seg1_max_Dev3 - Axis2_seg1_min_Dev3 # We expect that from a point and then, the length of the segment will change
                            #Be careful, Axis_2_face_Diff: this variable remains stable during the loop, it doesn't change
                            Face_seg_diff = abs(Axis_2_seg_Diff - Axis_2_face_Diff)
                            #Make the Algorithm stop some time.
                            if Step_sum >= Max_deviation:
                                Face_seg_diff = 0
                                Axis3_Edge_Z1_max_new = 0
                    Deviation_in_Max_in_axis3 = Axis3_Edge_Z1_max_new - Axis3_Edge_Z1_max_initial  #Deviation is extended in negative corresponding axis, because Axis3_Edge_Z1_max < Axis3_Edge_Z1_max_initial then always Deviation_in_Max_in_axis3 < 0, as we desire
                    face_colour = Colour_allocation_OE(Deviation_in_Max_in_axis3)
                    New_Colour_face_object[2] = face_colour
                    print("Deviation Axis 3 max checked")
                    

                #For bounding Box 2.2:
                if points_Bound_Box_2_2.size > 0: #this means that array is not empty, then it might have inward or outward deviation 
                   # The Bounding Box 2.2: is extended in Axis 2 between Axis2_Edge2_min and Axis2_Edge2_max,
                    #and in Axis 3 between Axis3_Edge_Z2_min and Axis3_Edge_Z2_max
                    Axis1_seg2_max, Axis2_seg2_max, Axis3_seg2_max, Axis1_seg2_min, Axis2_seg2_min, Axis3_seg2_min = Segemnt_Coordinates_in_BB(points_Bound_Box_2_2, face_OE)
                    #Now we have all the values that represent the points of the segment that are included inside each bounding box 2.2 Apparently we are not going to use all values.
                    #We want to maintain the values of these variables that is why we have to be very careful for their names to not be replaced
                    
                    #1) FOR AXIS 2
                    #For Deviations in Axis 2, and from the edge and inwards: 
                    #We find the lengths of Face in axis 3, and of the segment in Axis 3 for the points included in BB 2.2
                    Axis_3_seg_Diff = Axis3_seg2_max - Axis3_seg2_min #From initial points inside the BB 2.2
                    Axis_3_face_Diff = Axis3_Edge_Z2_max - Axis3_Edge_Z2_min
                    Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                    Deviation_in_Max_in_axis2 = 0
                    Axis2_Edge2_max_new = Axis2_Edge2_max
                    Axis2_Edge2_max_initial = Axis2_Edge2_max
                    points_Bound_Box_2_2_new = points_Bound_Box_2_2.copy()
                    if Face_seg_diff <= 0.05: #This means that the points are extended in the entire length of the face in BB 2.2 in Axis 2 so, we are sure that we have deviation in Z axis, as the entire wall has moved as a whole
                        Step_sum = 0
                        while Face_seg_diff <= 0.05:# We are looking for inward deviation
                            #in this loop we deduct0.005 m, namely 1 cm every time from the right boundary of the BB 2.2. With this way we deduct the size of the BB until the above condition not satisfied
                            #We are in BB2.2 so we substract 0.005 in the variable: Axis2_Edge2_max_new
                            Axis2_Edge2_max_new = Axis2_Edge2_max_new - 0.005
                            #Increase the step by 0.005 to make the algorithm stop if there is no points at all
                            Step_sum = Step_sum + 0.005
                            #then we recalculate the Bounding box by incorporating this adjustment
                            points_Bound_Box_2_2_new = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge2_min, Axis2_Edge2_max_new, Axis3_Edge_Z2_min, Axis3_Edge_Z2_max)
                            Axis1_seg2_max_Dev2, Axis2_seg2_max_Dev2, Axis3_seg2_max_Dev2, Axis1_seg2_min_Dev2, Axis2_seg2_min_Dev2, Axis3_seg2_min_Dev2 = Segemnt_Coordinates_in_BB(points_Bound_Box_2_2_new, face_OE)
                            
                            Axis_3_seg_Diff = Axis3_seg2_max_Dev2 - Axis3_seg2_min_Dev2 # We expect that from a point and then, the length of the segment will change
                            #Be careful, Axis_3_face_Diff remains stable during the loop, it doesn't change
                            Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                            #Make the Algorithm stop some time.
                            if Step_sum >= Max_deviation:
                                Face_seg_diff = 10
                                Axis2_Edge2_max_new = 0
                    else:
                        Step_sum = 0
                        while Face_seg_diff >= 0.05: # We are looking for outward deviation
                            #in this loop we add 0.005 m, namely 1 cm every time in the right boundary of the BB 2.2. With this way we extend the size of the BB until the above condition not satisfied
                            #We are in BB2.2 so we add 0.005 in the variable: Axis2_Edge2_max_new
                            Axis2_Edge2_max_new = Axis2_Edge2_max_new + 0.005
                            #Increase the step by 0.005 to make the algorithm stop if there is no points at all
                            Step_sum = Step_sum + 0.005
                            #then we recalculate the Bounding box by incorporating this adjustment
                            points_Bound_Box_2_2_new = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge2_min, Axis2_Edge2_max_new, Axis3_Edge_Z2_min, Axis3_Edge_Z2_max)
                            Axis1_seg2_max_Dev2, Axis2_seg2_max_Dev2, Axis3_seg2_max_Dev2, Axis1_seg2_min_Dev2, Axis2_seg2_min_Dev2, Axis3_seg2_min_Dev2 = Segemnt_Coordinates_in_BB(points_Bound_Box_2_2_new, face_OE)
                            
                            Axis_3_seg_Diff = Axis3_seg2_max_Dev2 - Axis3_seg2_min_Dev2 # We expect that from a point and then, the length of the segment will change
                            #Be careful, Axis_3_face_Diff remains stable during the loop, it doesn't change
                            Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                            #Make the Algorithm stop some time.
                            if Step_sum >= Max_deviation:
                                Face_seg_diff = 0
                                Axis2_Edge2_max_new = 0
                    Deviation_in_Max_in_axis2 = Axis2_Edge2_max_new - Axis2_Edge2_max_initial  #Deviation is extended in positive corresponding axis, because always Axis2_Edge2_max > Axis2_Edge2_max_initial, always we have Deviation_in_Min_in_axis2 > 0, as we desire 
                    face_colour = Colour_allocation_OE(Deviation_in_Max_in_axis2)
                    New_Colour_face_object[3] = face_colour
                    print("Deviation Axis 2 max checked")
                    #2) FOR AXIS 3 
                    print("Deviation Axis 3 min cannot be checked")
                else: #This means that the array is empty so there are no inward Deviations, in this case we check for outward deviations, if there is none then the OE edge is considered correct
                    Axis1_seg2_max, Axis2_seg2_max, Axis3_seg2_max, Axis1_seg2_min, Axis2_seg2_min, Axis3_seg2_min = Segemnt_Coordinates_in_BB(points_Bound_Box_2_2, face_OE)
                    #1) FOR AXIS 2
                    #For Deviations in Axis 2, and from the edge and inwards: 
                    #We find the lengths of Face in axis 3, and of the segment in Axis 3 for the points included in BB 2.2
                    Axis_3_seg_Diff = Axis3_seg2_max - Axis3_seg2_min #From initial points inside the BB 2.2
                    Axis_3_face_Diff = Axis3_Edge_Z2_max - Axis3_Edge_Z2_min
                    Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                    Deviation_in_Max_in_axis2 = 0
                    Axis2_Edge2_max_new = Axis2_Edge2_max
                    Axis2_Edge2_max_initial = Axis2_Edge2_max
                    points_Bound_Box_2_2_new = points_Bound_Box_2_2.copy()
                    if Face_seg_diff >= 0.05: #This means that the points are not extended in the entire length of the face in BB 2.2 in Axis 2 so, we are not sure what kind of deviation we have, we apply Bounding box regression until this condition is not satisfied, then we check the deviation that we have found if it is above 5cm then it is considerable deviation
                        Step_sum = 0
                        while Face_seg_diff >= 0.05: # We are looking for outward deviation
                            #in this loop we add 0.005 m, namely 1 cm every time in the right boundary of the BB 2.2. With this way we extend the size of the BB until the above condition not satisfied
                            #We are in BB2.2 so we add 0.005 in the variable: Axis2_Edge2_max_new
                            Axis2_Edge2_max_new = Axis2_Edge2_max_new + 0.005
                            #Increase the step by 0.005 to make the algorithm stop if there is no points at all
                            Step_sum = Step_sum + 0.005
                            #then we recalculate the Bounding box by incorporating this adjustment
                            points_Bound_Box_2_2_new = Bounding_box_and_point_selection(face_OE, points, Axis1_minus, Axis1_max, Axis2_Edge2_min, Axis2_Edge2_max_new, Axis3_Edge_Z2_min, Axis3_Edge_Z2_max)
                            Axis1_seg2_max_Dev2, Axis2_seg2_max_Dev2, Axis3_seg2_max_Dev2, Axis1_seg2_min_Dev2, Axis2_seg2_min_Dev2, Axis3_seg2_min_Dev2 = Segemnt_Coordinates_in_BB(points_Bound_Box_2_2_new, face_OE)
                            
                            Axis_3_seg_Diff = Axis3_seg2_max_Dev2 - Axis3_seg2_min_Dev2 # We expect that from a point and then, the length of the segment will change
                            #Be careful, Axis_3_face_Diff remains stable during the loop, it doesn't change
                            Face_seg_diff = abs(Axis_3_seg_Diff - Axis_3_face_Diff)
                            #Make the Algorithm stop some time.
                            if Step_sum >= Max_deviation:
                                Face_seg_diff = 0
                                Axis2_Edge2_max_new = 0
                    Deviation_in_Max_in_axis2 = Axis2_Edge2_max_new - Axis2_Edge2_max_initial  #Deviation is extended in positive corresponding axis, because always Axis2_Edge2_max > Axis2_Edge2_max_initial, always we have Deviation_in_Min_in_axis2 > 0, as we desire 
                    face_colour = Colour_allocation_OE(Deviation_in_Max_in_axis2)
                    New_Colour_face_object[3] = face_colour
                    print("Deviation Axis 2 max checked")
                    #2) FOR AXIS 3 
                    print("Deviation Axis 3 min cannot be checked")
                #We check if the array that contains the Colours derived from the different deviation checking, has any values in it. If it has it means that BB1 or BB2, or BB1 and BB2, include point cloud points, and therefore the final colour of the face will be defined by this array.
                # Check if all elements are None
                if np.any([x is not None for x in New_Colour_face_object]): #This means that the array with the colours is not empty and the Final colour of the Face will be identified by its values
                    print("At least one element is not None")
                    face_colour_final = [0,153,0] # We set the final face colour as green by default. Then we are looking if there is any red colour produced from the five types of deviation
                    for j in range(len(New_Colour_face_object)):
                        #if at least one of the colours are red then at least one within the five types of deviations is exceed the tolerance
                        if New_Colour_face_object[j] == [153,0,0]:
                            face_colour_final = New_Colour_face_object[j] #then the face is coloured Red.
                    #then write this face as exisitning face with its corresposning deviations
                    New_Final_face_object = np.array([New_FaceOE[0], New_FaceOE[1],New_FaceOE[2], New_FaceOE[3], New_FaceOE[4], New_FaceOE[5], New_FaceOE[6], New_FaceOE[7], New_FaceOE[8], New_FaceOE[9], New_FaceOE[10], New_FaceOE[11], New_FaceOE[12], New_FaceOE[13], New_FaceOE[14], New_FaceOE[15], New_FaceOE[16], New_FaceOE[17], New_FaceOE[18], New_FaceOE[19], New_FaceOE[20], New_FaceOE[21], New_FaceOE[22], New_FaceOE[23], New_FaceOE[24], New_FaceOE[25], New_FaceOE[26], New_FaceOE[27], New_FaceOE[28], face_colour_final, None, Deviation_in_Max_in_axis2, Deviation_in_Min_in_axis2, Deviation_in_Max_in_axis3, None ], dtype=object)  #Deviation_in_1_axis, Deviation_in_Max_in_axis2, Deviation_in_Min_in_axis2, Deviation_in_Max_in_axis3, Deviation_in_Min_in_axis3
                    Final_Face_Segment = np.vstack([Final_Face_Segment, New_Final_face_object])
                elif np.all([x is None for x in New_Colour_face_object]): #This means that the array is empty and that there are no deviations
                    print("All elements are None")
                    face_colour_final = face_colour #so the Face is coloured Green because it is concluded that there are no points inside the two boxes
                    New_Final_face_object = np.array([New_FaceOE[0], New_FaceOE[1],New_FaceOE[2], New_FaceOE[3], New_FaceOE[4], New_FaceOE[5], New_FaceOE[6], New_FaceOE[7], New_FaceOE[8], New_FaceOE[9], New_FaceOE[10], New_FaceOE[11], New_FaceOE[12], New_FaceOE[13], New_FaceOE[14], New_FaceOE[15], New_FaceOE[16], New_FaceOE[17], New_FaceOE[18], New_FaceOE[19], New_FaceOE[20], New_FaceOE[21], New_FaceOE[22], New_FaceOE[23], New_FaceOE[24], New_FaceOE[25], New_FaceOE[26], New_FaceOE[27], New_FaceOE[28], face_colour_final, None, Deviation_in_Max_in_axis2, Deviation_in_Min_in_axis2, Deviation_in_Max_in_axis3, None ], dtype=object)  #Deviation_in_1_axis, Deviation_in_Max_in_axis2, Deviation_in_Min_in_axis2, Deviation_in_Max_in_axis3, Deviation_in_Min_in_axis3
                    Final_Face_Segment = np.vstack([Final_Face_Segment, New_Final_face_object])
        if i != 2: 
            #See if the OE Face exist at least in the another PLANE 
            # Filter rows where "IfcType" is "IfcWall"
            filtered_rows_2nd = Corrrect_Mat_Seg_2nd_axis[1:][Corrrect_Mat_Seg_2nd_axis[1:, 2] == "IfcWall"]
            # Combine the header with the filtered rows
            IfcWall_rows_2nd = np.vstack((Corrrect_Mat_Seg_2nd_axis[0], filtered_rows_2nd))
            #We initialise an index that represents true or false
            boolean_index = 0
            for face_Wall_2nd in IfcWall_rows_2nd[1:]:
                if  face_Wall_2nd[5] < face_OE[5] < face_Wall_2nd[5] + 1: #Here we find the Wall face that is assigned to the current OpenElement face by checking if the OpeningElement's Local ID is based on Wall's face Local ID, Ex: Wallface ID = 1, then all its opening Elements should have IDs, 1.1, 1.2, 1.3, etc
                    face_colour = [0,153,0] #(Red, Green, Blue)
                    boolean_index = int(1)
                else:
                    print("Array contains only header row and the algorithm haven't identified any matching segment for this array")
                    print("The segment is missing so the OE face is missing")
                    face_colour = [153,0,0] #(Red, Green, Blue)
            if boolean_index == 1:
                OE_face_index = 1 #This means that at least the OE face matches with a wall but it is excluded from the analysis
            else:
                OE_face_index = 0 # this means that the OE face does not match with any wall and does not exist at all.
            New_Missing_face_object = np.array([face_OE[0], face_OE[1],face_OE[2], face_OE[3], face_OE[4], face_OE[5], face_OE[6], face_OE[7], face_OE[8], face_OE[9], face_OE[10], face_OE[11], face_OE[12], face_OE[13], face_OE[14], face_OE[15], face_OE[16], face_OE[17], face_colour, OE_face_index], dtype=object)
            Missing_Faces = np.vstack([Missing_Faces, New_Missing_face_object])

    #the segment is mising therefore the Wall is missing and therefore, the Opening element is missing, so add that to
    return Final_Face_Segment, Missing_Faces
input_Segment_directory_2 = output_Segment_directory #this Directory has came from the Segment analysis part
#For XZ plane in terms of GCS
Final_Face_SegmentsXZ_OE = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Xmax", "Xmin", "Zmax", "Zmin", "Colour", "Deviation in Y-axis", "Deviation of face Xmax", "Deviation of face Xmin",  "Deviation of face Zmax", "Deviation of face Zmin" ]], dtype=object) # the axes in these arrays refer to the GCS
Missing_Faces_XZ_OE = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin", "Colour", "Label_index" ]], dtype=object)
Final_Face_SegmentsXZ_OE, Missing_Faces_XZ_OE = OpeningElement_Comparison(Correct_Matching_SegmentsXZ, Correct_Matching_SegmentsYZ, FacesXZ_OE, Missing_Faces_XZ_OE, Final_Face_SegmentsXZ_OE, input_Segment_directory_2)
print(Final_Face_SegmentsXZ_OE)
#For YZ plane in terms of GCS
Final_Face_SegmentsYZ_OE = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Ymax", "Ymin", "Zmax", "Zmin", "Colour", "Deviation in X-axis", "Deviation of face Ymax", "Deviation of face Ymin",  "Deviation of face Zmax", "Deviation of face Zmin" ]], dtype=object) #35 Columns, 34 indices
Missing_Faces_YZ_OE = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin", "Colour", "Label_index" ]], dtype=object)
Final_Face_SegmentsYZ_OE, Missing_Faces_YZ_OE = OpeningElement_Comparison(Correct_Matching_SegmentsYZ, Correct_Matching_SegmentsXZ,  FacesYZ_OE, Missing_Faces_YZ_OE, Final_Face_SegmentsYZ_OE, input_Segment_directory_2)
print(Final_Face_SegmentsYZ_OE)
#Merge all the results under one array for every plane XZ, YZ, XY, based on the structure of the arrays FacesXZ, FacesYZ, FacesXY
def Merging_the_result_arrays(Faces, Excluded_Faces, Missing_Faces, Final_Face_Segments, Merged_faces):
    for face in Faces[1:]:
        Missing_Label = "Missing Face"
        Exisitng_Label = "Exists"
        Excluded_label = "Excluded"
        for excl_face in Excluded_Faces[1:]:
            if face[1] == excl_face[1] and face[2] == excl_face[2] and face[3] == excl_face[3]:
                face_colour = [255,255,102] #Yellow becasue they are not included in the analysis, becasue there is no Segment to be assigned to them 
                New_Merged_face_object = np.array([excl_face[0], excl_face[1],excl_face[2], excl_face[3], excl_face[4], excl_face[5], excl_face[6], excl_face[7], excl_face[8], excl_face[9], excl_face[10], excl_face[11], excl_face[12], excl_face[13], excl_face[14], excl_face[15], excl_face[16], excl_face[17], None, None, None, None, None, None, None, None, None, None, None, face_colour, None, None, None, None, None, excl_face[13], Excluded_label], dtype=object)
                Merged_faces = np.vstack([Merged_faces, New_Merged_face_object])
        for missing_face in Missing_Faces[1:]:
            if face[1] == missing_face[1] and face[2] == missing_face[2] and face[3] == missing_face[3]:
                face_colour = [1,0,0] #Red becasue they are missing from the model, the haven't been built yet
                New_Merged_face_object = np.array([missing_face[0], missing_face[1],missing_face[2], missing_face[3], missing_face[4], missing_face[5], missing_face[6], missing_face[7], missing_face[8], missing_face[9], missing_face[10], missing_face[11], missing_face[12], missing_face[13], missing_face[14], missing_face[15], missing_face[16], missing_face[17], None, None, None, None, None, None, None, None, None, None, None, missing_face[18], None, None, None, None, None, missing_face[13], Missing_Label ], dtype=object)
                Merged_faces = np.vstack([Merged_faces, New_Merged_face_object])
        for Final_face in Final_Face_Segments[1:]:
            if face[1] == Final_face[1] and face[2] == Final_face[2] and face[3] == Final_face[3]:
                face_colour = [0,1,0] #Green becasue the are not included in the analysis and they are deemed correct
                New_Merged_face_object = np.array([Final_face[0], Final_face[1],Final_face[2], Final_face[3], Final_face[4], Final_face[5], Final_face[6], Final_face[7], Final_face[8], Final_face[9], Final_face[10], Final_face[11], Final_face[12], Final_face[13], Final_face[14], Final_face[15], Final_face[16], Final_face[17], Final_face[18], Final_face[19], Final_face[20], Final_face[21], Final_face[22], Final_face[23], Final_face[24], Final_face[25], Final_face[26], Final_face[27], Final_face[28], Final_face[29], Final_face[30], Final_face[31], Final_face[32], Final_face[33], Final_face[34], Final_face[13], Exisitng_Label], dtype=object)
                Merged_faces = np.vstack([Merged_faces, New_Merged_face_object])
    return Merged_faces
Merged_facesXZ = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: y", "Dimension: dx", "Dimension: dz", "Average Location: x", "Average Location: z", "Xmax", "Xmin", "Zmax", "Zmin", "Colour", "Deviation in Y-axis", "Deviation of face Xmax", "Deviation of face Xmin",  "Deviation of face Zmax", "Deviation of face Zmin", "Adjusted vertices", "Existance label" ]], dtype=object) #The new Merged array has 34 columns, for safety we include all the data that we have calculated during the analysis in this Final array
Merged_facesXZ = Merging_the_result_arrays(FacesXZ, Excluded_FaceXZ, Missing_FacesXZ, Final_Face_SegmentXZ, Merged_facesXZ)
        
Merged_facesYZ = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Face Normal", "Face Vertices", "Ymax", "Ymin", "Zmax", "Zmin", "Segment ID", "Segment Name", "Segment Location: x", "Dimension: dy", "Dimension: dz", "Average Location: y", "Average Location: z", "Ymax", "Ymin", "Zmax", "Zmin", "Colour", "Deviation in X-axis", "Deviation of face Ymax", "Deviation of face Ymin",  "Deviation of face Zmax", "Deviation of face Zmin", "Adjusted vertices", "Existance label" ]], dtype=object) #35 Columns, 34 indices
Merged_facesYZ = Merging_the_result_arrays(FacesYZ, Excluded_FaceYZ, Missing_FacesYZ, Final_Face_SegmentYZ, Merged_facesYZ)    

Merged_facesXY = np.array([["Face Type in GCS", "Face Type in LCS", "IfcType", "GUID", "ID", "Local ID", "Floor Level (z)","Face Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Face Normal", "Face Vertices", "Xmax", "Xmin", "Ymax", "Ymin", "Segment ID", "Segment Name", "Segment Location: z", "Dimension: dx", "Dimension: dy", "Average Location: x", "Average Location: y", "Xmax", "Xmin", "Ymax", "Ymin", "Colour", "Deviation in Z-axis", "Deviation of face Xmax", "Deviation of face Xmin",  "Deviation of face Ymax", "Deviation of face Ymin", "Adjusted vertices", "Existance label" ]], dtype=object)
Merged_facesXY = Merging_the_result_arrays(FacesXY, Excluded_FaceXY, Missing_FacesXY, Final_Face_SegmentXY, Merged_facesXY) 
print(Merged_facesXZ)
def Adjust_face_vertices(face_OE): #This works only for visualisation purposes in order the objects to be more clear and not folded one each other
    if face_OE[0]=="Face XZ": #If the face is located in XZ plane the Axis 1, 2, 3 is the Y, X, Z axes respectively
        #Then the common coordinate of the all the face vertices is the y.
        # Extract the sign of the first coordinate
        Face_Normal = face_OE[12]
        if Face_Normal[1] < 0:
            first_coordinate_sign = -1
        elif Face_Normal[1] > 0:
            first_coordinate_sign = +1
        adjusted_vertices = face_OE[13].copy()
        adjustment = 0.02*first_coordinate_sign
        # Add the adjustment to the x coordinates
        for i in range(len(adjusted_vertices)):
            adjusted_vertices[i, 1] =  adjusted_vertices[i, 1] + adjustment
    elif face_OE[0] =="Face YZ":
        #Then the common coordinate of the all the face vertices is the x.
        # Extract the sign of the first coordinate
        Face_Normal = face_OE[12]
        if Face_Normal[0] < 0:
            second_coordinate_sign = -1
        elif Face_Normal[0] > 0:
            second_coordinate_sign = +1
        adjusted_vertices = face_OE[13].copy()
        adjustment = 0.02*second_coordinate_sign
        # Add the adjustment to the x coordinates
        for i in range(len(adjusted_vertices)):
            adjusted_vertices[i, 0] =  adjusted_vertices[i, 0] + adjustment
    return adjusted_vertices
def Merging_the_result_arrays_OE(Faces_OE, Missing_Faces_OE, Final_Faces_OE, Merged_faces):
    Missing_Label = "Missing Face"
    Exisitng_Label = "Exists"
    Excluded_label = "Excluded"
    for face in Faces_OE[1:]:
        for missing_face_OE in Missing_Faces_OE[1:]:
            if face[1] == missing_face_OE[1] and face[2] == missing_face_OE[2] and face[3] == missing_face_OE[3]:
                adjusted_vertices = Adjust_face_vertices(missing_face_OE)
                if missing_face_OE[19] == int(1): 
                    OE_label = Excluded_label
                elif missing_face_OE[19] == int(0): 
                    OE_label = Missing_Label
                New_Merged_face_object = np.array([missing_face_OE[0], missing_face_OE[1],missing_face_OE[2], missing_face_OE[3], missing_face_OE[4], missing_face_OE[5], missing_face_OE[6], missing_face_OE[7], missing_face_OE[8], missing_face_OE[9], missing_face_OE[10], missing_face_OE[11], missing_face_OE[12], missing_face_OE[13], missing_face_OE[14], missing_face_OE[15], missing_face_OE[16], missing_face_OE[17], None, None, None, None, None, None, None, None, None, None, None, missing_face_OE[18], None, None, None, None, None, adjusted_vertices, OE_label ], dtype=object)
                Merged_faces = np.vstack([Merged_faces, New_Merged_face_object])
        for Final_Face_OE in Final_Faces_OE[1:]:
            if face[1] == Final_Face_OE[1] and face[2] == Final_Face_OE[2] and face[3] == Final_Face_OE[3]:
                adjusted_vertices = Adjust_face_vertices(Final_Face_OE)
                New_Merged_face_object =  np.array([Final_Face_OE[0], Final_Face_OE[1],Final_Face_OE[2], Final_Face_OE[3], Final_Face_OE[4], Final_Face_OE[5], Final_Face_OE[6], Final_Face_OE[7], Final_Face_OE[8], Final_Face_OE[9], Final_Face_OE[10], Final_Face_OE[11], Final_Face_OE[12], Final_Face_OE[13], Final_Face_OE[14], Final_Face_OE[15], Final_Face_OE[16], Final_Face_OE[17], Final_Face_OE[18], Final_Face_OE[19], Final_Face_OE[20], Final_Face_OE[21], Final_Face_OE[22], Final_Face_OE[23], Final_Face_OE[24], Final_Face_OE[25], Final_Face_OE[26], Final_Face_OE[27], Final_Face_OE[28], Final_Face_OE[29], Final_Face_OE[30], Final_Face_OE[31], Final_Face_OE[32], Final_Face_OE[33], Final_Face_OE[34], adjusted_vertices, Exisitng_Label ], dtype=object)
                Merged_faces = np.vstack([Merged_faces, New_Merged_face_object])
    return Merged_faces

Merged_facesXZ = Merging_the_result_arrays_OE(FacesXZ_OE, Missing_Faces_XZ_OE, Final_Face_SegmentsXZ_OE, Merged_facesXZ)        
Merged_facesYZ = Merging_the_result_arrays_OE(FacesYZ_OE, Missing_Faces_YZ_OE, Final_Face_SegmentsYZ_OE, Merged_facesYZ)    
#Creating OBJ files for Visualising the results in MeshLab
def create_a_mesh(face, individual_triangle):
    # Create a mesh object
    mesh = o3d.geometry.TriangleMesh()

    # Assuming face is a list of vertices of the current face
    mesh.vertices = o3d.utility.Vector3dVector(face)
    mesh.triangles = o3d.utility.Vector3iVector(individual_triangle)
    return(mesh)
def create_a_mesh_with_colour(face, individual_triangle, color_face_indices, face_colour):
    # Create a mesh object
    mesh = o3d.geometry.TriangleMesh()

    # Assuming face is a list of vertices of the current face
    mesh.vertices = o3d.utility.Vector3dVector(face)
    mesh.triangles = o3d.utility.Vector3iVector(individual_triangle)
    # If color_face_indices is provided, color those faces with red
    if color_face_indices is not None:
        num_vertices = len(face) # OBJ file format assign colours to individual vertices, so we need the vertices that constitute the face
        colors = np.zeros((num_vertices, 3))  # Initialize colors array, this is a 4x3 array, inside which we store the colour array of its vertex

        # Assign the predifined colour of the face, "face_colour" variable in every vertex in the face, first position=colour of vertex 1, second position for v2,etc
        for i in color_face_indices:
            colors[i] = face_colour

        mesh.vertex_colors = o3d.utility.Vector3dVector(colors)
    return(mesh)
individual_triangle = np.array([
    [0, 1, 2 ],  # Tr1 
    [2, 3, 0 ],  # Tr2 
    ])
Combine_meshes = []
# Iterating through the array
for row in Merged_facesXZ[1:]:
    # Iterate over each face and create a mesh
    if row[3] == Higher_slab_GUId:
        print("This slab is excluded from the visualisation, because constitute the upper slab of the floor")
    else:
        face = row[35]  #row[13] for real face vertices
        face_colour = row[29]
        # Create a mesh object
        mesh = create_a_mesh_with_colour(face, individual_triangle, individual_triangle, face_colour)
        # Append the mesh to the list
        Combine_meshes.append(mesh) #For creating the Combine mesh
for row in Merged_facesYZ[1:]:
    # Iterate over each face and create a mesh
    if row[3] == Higher_slab_GUId:
        print("This slab is excluded from the visualisation, because constitute the upper slab of the floor")
    else:
        face = row[35]  #row[13] for real face vertices
        face_colour = row[29]
        # Create a mesh object
        mesh = create_a_mesh_with_colour(face, individual_triangle, individual_triangle, face_colour)
        # Append the mesh to the list
        Combine_meshes.append(mesh) #For creating the Combine mesh
for row in Merged_facesXY[1:]:
    # Iterate over each face and create a mesh
    if row[3] == Higher_slab_GUId:
        print("This slab is excluded from the visualisation, because constitute the upper slab of the floor")
    else:
        face = row[35]  #row[13] for real face vertices
        face_colour = row[29]
        # Create a mesh object
        mesh = create_a_mesh_with_colour(face, individual_triangle, individual_triangle, face_colour)
        # Append the mesh to the list
        Combine_meshes.append(mesh) #For creating the Combine mesh

def merge_meshes_into_one_mesh(ifcObject_mesh): # This function merge the seperate meshes into one mesh in order to be visualised properly in Meshlab tool
    merged_mesh = o3d.geometry.TriangleMesh()
    for mesh in ifcObject_mesh:
        merged_mesh += mesh
    return merged_mesh
Combine_Merged_meshes = merge_meshes_into_one_mesh(Combine_meshes)
Combine_merged_mesh_file_path = r"C:/Python Output/Real data/OBJ files/Final/Combined_Geometry_SDK_real_Final_7.obj"
# Export combined geometry to OBJ or PLY
o3d.io.write_triangle_mesh(Combine_merged_mesh_file_path, Combine_Merged_meshes)

def convert_to_csv(faces_array, csv_file_path):
    # Convert numpy array to pandas DataFrame
    df = pd.DataFrame(faces_array[1:], columns=faces_array[0])
    # Write to CSV
    df.to_csv(csv_file_path, index=False)

csv_file_pathXZ = r"C:/Python Output/Real data/CSV files/Final/Final_XZ_7.csv"
convert_to_csv(Merged_facesXZ, csv_file_pathXZ)
csv_file_pathYZ = r"C:/Python Output/Real data/CSV files/Final/Final_YZ_7.csv"
convert_to_csv(Merged_facesYZ, csv_file_pathYZ)
csv_file_pathXY = r"C:/Python Output/Real data/CSV files/Final/Final_XY_7.csv"
convert_to_csv(Merged_facesXY, csv_file_pathXY)