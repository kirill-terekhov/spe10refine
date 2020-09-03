#include "inmost.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <amesh.h>

using namespace INMOST;


double * perm = NULL; 
double * por = NULL;
std::string problem_name = "spe10_inclined";

#define ID(i,j,k) ((i) + (j)*60 + (k)*220*60)
#define V_ID(x, y, z) ((x-x_beg)*(y_end-y_beg+1)*(z_end-z_beg+1) + (y-y_beg)*(z_end-z_beg+1) + (z-z_beg))


class CopySpe10Data : public AbstractSubModel
{
	TagReal phi;
	TagRealArray perm;
public:
	CopySpe10Data(Tag phi, Tag perm) : phi(phi), perm(perm) {}
	bool PrepareEntries(Model& m) {return true;}
	bool Initialize(Model& m) {return true;}
	bool FillResidual(Residual& R) const {return true;}
	bool UpdateSolution(const Sparse::Vector& V, double alpha) {return true;}
	bool UpdateTimeStep() {return true;}
	bool SetTimeStep(double dt) {return true;}
	bool RestoreTimeStep() {return true;}
	void Adaptation(Mesh & m, SearchKDTree & tree) const 
	{
		rMatrix xn(3,1);
		for(Mesh::iteratorCell it = m.BeginCell(); it != m.EndCell(); ++it)
		{
			Cell c = it->self();
			if( c.New() )
			{
				c.Centroid(xn.data());
				Cell old = tree.SearchCell(xn.data());
				if( old.isValid() )
				{
					phi[c] = phi[old];
					perm(c,3,1) = perm(old,3,1);
				}
				else
				{
					std::cout << "cannot find old cell for new cell " << it->LocalID() << std::endl;
					old = tree.SearchCell(xn.data(),true);
				}
			}
		}
	}
};

bool read_spe10()
{
	perm = new double[60*220*85*3];
	por = new double[60*220*85];
	if( por == NULL || perm == NULL )
	{
		std::cout << "Not enough memory for spe10 data." << std::endl;
		return false;
	}
	FILE * fperm = fopen("spe_perm.dat","r");
	FILE * fpor = fopen("spe_phi.dat","r");
	if( fperm == NULL || fpor == NULL ) 
	{
		std::cout << "SPE10 data files not found." << std::endl;
		std::cout << "Expecting spe_perm.dat and spe_phi.dat." << std::endl;
		return false;
	}

	for(int k = 0; k < 85; k++)
	{
		for(int j = 0; j < 220; j++)
		{
			for(int i = 0; i < 60; i++)
			{
				if(fscanf(fpor,"%lg ", &por[ID(i,j,k)]) != 1 )
				{
					std::cout << "error in phi!" << std::endl;
					return false;
				}
			}
		}
	}
	for(int k = 0; k < 85; k++)
	{
		for(int j = 0; j < 220; j++)
		{
			for(int i = 0; i < 60; i++)
			{
				if(fscanf(fperm,"%lg ",&perm[ID(i,j,k)*3+0]) != 1 )
				{
					std::cout << "error in perm!" << std::endl;
					return false;
				}
			}
		}
	}
	for(int k = 0; k < 85; k++)
	{
		for(int j = 0; j < 220; j++)
		{
			for(int i = 0; i < 60; i++)
			{
				if(fscanf(fperm,"%lg ",&perm[ID(i,j,k)*3+1]) != 1 )
				{
					std::cout << "error in perm!" << std::endl;
					return false;
				}
			}
		}
	}

	for(int k = 0; k < 85; k++)
	{
		for(int j = 0; j < 220; j++)
		{
			for(int i = 0; i < 60; i++)
			{
				if(fscanf(fperm,"%lg ",&perm[ID(i,j,k)*3+2]) != 1 )
				{
					std::cout << "error in perm!" << std::endl;
					return false;
				}
			}
		}
	}


	fclose(fperm);
	fclose(fpor);
	return true;
}

void kill_spe10()
{
	if( perm != NULL ) delete [] perm;
	if( por != NULL ) delete [] por;
}


void CreateCubeElement(Mesh *m, ElementArray<Node> verts, ElementArray<Cell> & cells)
{
	// Define six cube faces assuming verts are numerated in the way presented above
	const INMOST_DATA_INTEGER_TYPE face_nodes[24] = {0,4,6,2, 1,3,7,5, 0,1,5,4, 2,6,7,3, 0,2,3,1, 4,5,7,6};
	const INMOST_DATA_INTEGER_TYPE num_nodes[6]   = {4,       4,       4,       4,       4,       4};

	cells.push_back(m->CreateCell(verts,face_nodes,num_nodes,6).first); // Create the cubic cell in the mesh
}

void CreateNWTetElements(Mesh *m, ElementArray<Node> verts, ElementArray<Cell> & cells)
{
	// Define prism faces assuming verts are numerated in the way presented above
	const INMOST_DATA_INTEGER_TYPE ne_face_nodes1[12] = {0,1,5,  5,1,3,   1,0,3, 3,0,5};
	const INMOST_DATA_INTEGER_TYPE ne_num_nodes1[4]   = {3,3,3,3};

	const INMOST_DATA_INTEGER_TYPE ne_face_nodes2[12] = {0,3,5, 0,7,3, 5,3,7, 0,5,7};
	const INMOST_DATA_INTEGER_TYPE ne_num_nodes2[4]   = {3,3,3,3};

	const INMOST_DATA_INTEGER_TYPE ne_face_nodes3[12] = {0,7,5, 4,5,7, 0,5,4, 0,4,7};
	const INMOST_DATA_INTEGER_TYPE ne_num_nodes3[4]   = {3,3,3,3};
	

	const INMOST_DATA_INTEGER_TYPE sw_face_nodes1[12] = {0,3,7, 2,7,3, 0,7,2, 0,2,3};
	const INMOST_DATA_INTEGER_TYPE sw_num_nodes1[4]   = {3,3,3,3};

	const INMOST_DATA_INTEGER_TYPE sw_face_nodes2[12] = {0,7,4, 0,2,7, 2,4,7, 0,4,2};
	const INMOST_DATA_INTEGER_TYPE sw_num_nodes2[4]   = {3,3,3,3};

	const INMOST_DATA_INTEGER_TYPE sw_face_nodes3[12] = {4,6,2, 6,7,2, 4,7,6, 4,2,7};
	const INMOST_DATA_INTEGER_TYPE sw_num_nodes3[4]   = {3,3,3,3};

	cells.push_back(m->CreateCell(verts,ne_face_nodes1,ne_num_nodes1,4).first); // Create north-east prismatic cell
	cells.push_back(m->CreateCell(verts,ne_face_nodes2,ne_num_nodes2,4).first); // Create north-east prismatic cell
	cells.push_back(m->CreateCell(verts,ne_face_nodes3,ne_num_nodes3,4).first); // Create north-east prismatic cell
	cells.push_back(m->CreateCell(verts,sw_face_nodes1,sw_num_nodes1,4).first); // Create south-west prismatic cell
	cells.push_back(m->CreateCell(verts,sw_face_nodes2,sw_num_nodes2,4).first); // Create south-west prismatic cell
	cells.push_back(m->CreateCell(verts,sw_face_nodes3,sw_num_nodes3,4).first); // Create south-west prismatic cell
}


int main(int argc, char *argv[]) 
{
	Mesh::Initialize(&argc,&argv);
#if defined(USE_PARTITIONER)
    Partitioner::Initialize(&argc,&argv); // Initialize the partitioner activity
#endif
	
	if (argc < 2)
	{
		std::cout << "Usage: " << argv[0] ;
		std::cout << " [incline=0] [tetrahedral=0] [refines=1] [x_beg=0] [x_end=60] [y_beg=0] [y_end=220] [z_beg=0] [z_end=85] [scale_permiability_x=1] [scale_permiability_y=1] [scale_permiability_z=1]" << std::endl;
		std::cout << " 0 <= x_beg < x_end <= 60 " << std::endl;
		std::cout << " 0 <= y_beg < y_end <= 220 " << std::endl;
		std::cout << " 0 <= z_beg < z_end <= 85 " << std::endl;
		return -1;
	}
	int tetra = 0;
	Mesh * mesh;
	double scale_permiability_x,scale_permiability_y,scale_permiability_z;
	int x_beg = 0, x_end = 60, y_beg = 0, y_end = 220, z_beg = 0, z_end = 85, refines = 1;
	bool incline = false;
	scale_permiability_x = scale_permiability_y=scale_permiability_z = 1;


	if( argc > 3 ) refines = atoi(argv[3]);
		

	mesh = new Mesh();
	mesh->SetCommunicator(INMOST_MPI_COMM_WORLD);
	
	Tag tporo = mesh->CreateTag("PORO",DATA_REAL,CELL,NONE,1);
	Tag tK = mesh->CreateTag("K",DATA_REAL,CELL,NONE,3);
	
	//only one processor builds the mesh
	if( mesh->GetProcessorRank() == 0 )
	{
		std::cout << "Read in spe10 data." << std::endl;
		if( !read_spe10() )
		{
			kill_spe10();
			return -1;
		}
		std::cout << "Done." << std::endl;

		

		if( argc > 1 ) incline = atoi(argv[1]);
		if( argc > 2) tetra = atoi(argv[2]);
		if( argc > 4 ) x_beg = atoi(argv[4]);
		if( argc > 5 ) x_end = atoi(argv[5]);
		if( argc > 6 ) y_beg = atoi(argv[6]);
		if( argc > 7 ) y_end = atoi(argv[7]);
		if( argc > 8 ) z_beg = atoi(argv[8]);
		if( argc > 9 ) z_end = atoi(argv[9]);
		if( argc > 10 ) scale_permiability_x = atof(argv[10]);
		if( argc > 11 ) scale_permiability_y = atof(argv[11]);
		if( argc > 12) scale_permiability_z = atof(argv[12]);
		
		if( tetra ) std::cout << "Building tetrahedral grid " << argv[2] << std::endl;
		else std::cout << "Building cubical grid " << argv[2] << std::endl;
		
		
		std::cout << "Create mesh nodes." << std::endl;


		for (int i = x_beg; i <= x_end; i++) 
		{
			for (int j = y_beg; j <= y_end; j++) 
			{
				for (int k = z_beg; k <= z_end; k++) 
				{
					Storage::real xyz[3];
					bool mark = false;
					xyz[0] = 240.0 * i / 60.0;
					xyz[1] = 440.0 * j / 220.0;
					xyz[2] = 340.0 * k / 85.0;

					if( incline )
					{
						double alpha = xyz[1]/440.0;
						if( xyz[0] > 120.0 )
							xyz[0] += (1.0-alpha)*(120-xyz[0]*0.5) + (alpha)*(80 + xyz[0]*0.5);
						else
							xyz[0] += (1.0-alpha)*(xyz[0]*0.5) + (alpha)*(200.0 - xyz[0]*0.5);
					}

					Node c = mesh->CreateNode(xyz);
					if (c->LocalID() != V_ID(i, j, k)) std::cout << "v_id = " << c->LocalID() << ", [i,j,k] = " << V_ID(i, j, k) << std::endl;
				}
			}
		}
		
		std::cout << "Done." << std::endl;
		
		std::cout << "Create mesh cells." << std::endl;

		Tag materials = mesh->CreateTag("MATERIAL",DATA_INTEGER,CELL,NONE,1);
		const INMOST_DATA_INTEGER_TYPE nvf[24] = { 2, 3, 1, 0, 4, 5, 7, 6, 0, 1, 5, 4, 3, 2, 6, 7, 2, 0, 4, 6, 1, 3, 7, 5 };
		const INMOST_DATA_INTEGER_TYPE numnodes[6] = { 4, 4, 4, 4, 4, 4 };
		for (int i = x_beg; i < x_end; i++) 
		{
			for (int j = y_beg; j < y_end; j++) 
			{
				for (int k = z_beg; k < z_end; k++) 
				{
					double poro = por[ID(i,j,k)];
					double Kabs = sqrt(perm[ID(i,j,k)*3+0]*perm[ID(i,j,k)*3+0]+perm[ID(i,j,k)*3+1]*perm[ID(i,j,k)*3+1]+perm[ID(i,j,k)*3+2]*perm[ID(i,j,k)*3+2]);
					if( poro > 1.0e-6 && Kabs > 1.0e-6)
					{
						/*
						if( perm[ID(i,j,k)*3+0] < 1.0e-6 ||
						perm[ID(i,j,k)*3+1] < 1.0e-6 ||
						perm[ID(i,j,k)*3+2] < 1.0e-6)
						{
						i,j,k,perm[ID(i,j,k)*3+0],perm[ID(i,j,k)*3+1],perm[ID(i,j,k)*3+2],poro);
						perm[ID(i,j,k)*3+0] = std::max(1.0e-6,perm[ID(i,j,k)*3+0]);
						perm[ID(i,j,k)*3+1] = std::max(1.0e-6,perm[ID(i,j,k)*3+1]);
						perm[ID(i,j,k)*3+2] = std::max(1.0e-6,perm[ID(i,j,k)*3+2]);
						}
						*/

						ElementArray<Node> verts(mesh,8); 
						verts[0] = mesh->NodeByLocalID(V_ID(i + 0, j + 0, k + 0));
						verts[1] = mesh->NodeByLocalID(V_ID(i + 1, j + 0, k + 0));
						verts[2] = mesh->NodeByLocalID(V_ID(i + 0, j + 1, k + 0));
						verts[3] = mesh->NodeByLocalID(V_ID(i + 1, j + 1, k + 0));
						verts[4] = mesh->NodeByLocalID(V_ID(i + 0, j + 0, k + 1));
						verts[5] = mesh->NodeByLocalID(V_ID(i + 1, j + 0, k + 1));
						verts[6] = mesh->NodeByLocalID(V_ID(i + 0, j + 1, k + 1));
						verts[7] = mesh->NodeByLocalID(V_ID(i + 1, j + 1, k + 1));
						ElementArray<Cell> cells;
						if( tetra )
							CreateNWTetElements(mesh,verts,cells);
						else
							CreateCubeElement(mesh,verts,cells);
						for(int q= 0; q < cells.size(); ++q)
						{
							Cell c = cells[q];
							c->RealDF(tporo) = poro;
							Storage::real_array tensor = c->RealArrayDF(tK);
							tensor[0] = perm[ID(i,j,k)*3+0]*scale_permiability_x;
							tensor[1] = perm[ID(i,j,k)*3+1]*scale_permiability_y;
							tensor[2] = perm[ID(i,j,k)*3+2]*scale_permiability_z;
						

							if( k > 50 ) c->IntegerDF(materials) = 1;
							else c->IntegerDF(materials) = 0;
						}
					}
				}
			}
		}
		
		std::cout << "Done." << std::endl;
		
		std::cout << "Delete empty cells." << std::endl;

		for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
		{
			ElementArray<Element> cells = it->BridgeAdjacencies(FACE,CELL);
			if( cells.empty() )
			{
				it->Delete();
			}
		}
		
		std::cout << "Done." << std::endl;
		
		std::cout << "Delete orphan elements." << std::endl;

		for(ElementType etype = FACE; etype >= NODE; etype = etype >> 1)
		{
			for(Mesh::iteratorElement it = mesh->BeginElement(etype); it != mesh->EndElement(); ++it)
				if( it->nbAdjElements(CELL) == 0 )
					it->Delete();
		}
		
		std::cout << "Done." << std::endl;

		{
			Storage::bulk_array name = mesh->self()->BulkArray(mesh->CreateTag("PROBLEMNAME",DATA_BULK,MESH,NONE));
			name.replace(name.begin(),name.end(),problem_name.begin(),problem_name.end());
		}


		std::cout << "Nodes: " << mesh->NumberOfNodes() << "." << std::endl;
		std::cout << "Edges: " << mesh->NumberOfEdges() << "." << std::endl;
		std::cout << "Faces: " << mesh->NumberOfFaces() << "." << std::endl;
		std::cout << "Cells: " << mesh->NumberOfCells() << "." << std::endl;

		std::cout << "I'm ready!" << std::endl;
	}
	//distribute the mesh
	
	
#if defined(USE_PARTITIONER)
	if( mesh->GetProcessorsNumber() > 1 )
	{
		if( !mesh->GetProcessorRank() ) std::cout << "Resolve shared elements." << std::endl;
		mesh->ResolveShared();
		if( !mesh->GetProcessorRank() ) std::cout << "Done." << std::endl;
		if( !mesh->GetProcessorRank() ) std::cout << "Partition the mesh." << std::endl;
		Partitioner p(mesh);
		p.SetMethod(Partitioner::INNER_KMEANS,Partitioner::Partition);
		//p.SetMethod(Partitioner::Parmetis,Partitioner::Partition);
		//p.SetMethod(Partitioner::Parmetis,Partitioner::Refine);
		p.Evaluate();
		if( !mesh->GetProcessorRank() ) std::cout << "Done." << std::endl;
		if( !mesh->GetProcessorRank() ) std::cout << "Redestribute the mesh." << std::endl;
		mesh->Redistribute();
		if( !mesh->GetProcessorRank() ) std::cout << "Done." << std::endl;
		if( !mesh->GetProcessorRank() ) std::cout << "Compactify data." << std::endl;
		mesh->ReorderEmpty(CELL|FACE|EDGE|NODE);
		if( !mesh->GetProcessorRank() ) std::cout << "Done." << std::endl;
		std::cout << "Number of cells on " << mesh->GetProcessorRank() << " is " << mesh->NumberOf(CELL) << std::endl;
		//mesh->Barrier();
	}
#endif
	//refine
	if( refines > 0 )
	{
		Automatizator aut;
		Model transdata(aut);
		CopySpe10Data copy(tporo,tK);
		transdata.AddSubModel("copy",copy);
		AdaptiveMesh refiner(*mesh);
		refiner.SetModel(&transdata);
		TagInteger indicator = mesh->CreateTag("indicator",DATA_INTEGER,CELL,NONE,1);
		for(int q = 0; q < refines; ++q)
		{
			if( !mesh->GetProcessorRank() ) std::cout << "Refinement cycle " << q+1 << "/" << refines << "." << std::endl;
			//all cells to be refined
			for(Mesh::iteratorCell it = mesh->BeginCell(); it != mesh->EndCell(); ++it)
				indicator[*it] = 1;
			if( !mesh->GetProcessorRank() ) std::cout << "Start refine." << std::endl;
			double t = Timer();
			refiner.Refine(indicator);
			t = Timer() - t;
			if( !mesh->GetProcessorRank() ) std::cout << "Done in " << t << "s." << std::endl;
			unsigned ncells = mesh->TotalNumberOf(CELL);
			unsigned nnodes = mesh->TotalNumberOf(NODE);
			unsigned nedges = mesh->TotalNumberOf(EDGE);
			unsigned nfaces = mesh->TotalNumberOf(FACE);
			if( !mesh->GetProcessorRank() ) 
			{
				std::cout << "Nodes: " << nnodes << "." << std::endl;
				std::cout << "Edges: " << nedges << "." << std::endl;
				std::cout << "Faces: " << nfaces << "." << std::endl;
				std::cout << "Cells: " << ncells << "." << std::endl;
			}
			//mesh->Barrier();
		}
	}
	if( !mesh->GetProcessorRank() ) std::cout << "Saving the grid into \"grid" << (mesh->GetProcessorsNumber() > 1 ? ".pvtk" : ".vtk") << "\"." << std::endl;
	if( mesh->GetProcessorsNumber() > 1 )
		mesh->Save("grid.pvtk");
	else mesh->Save("grid.vtk");
	if( !mesh->GetProcessorRank() ) std::cout << "File written!" << std::endl;

	kill_spe10();
	delete mesh;
	
#if defined(USE_PARTITIONER)
    Partitioner::Finalize(); // Finalize the partitioner activity
#endif
	Mesh::Finalize();
	return 0;
}
