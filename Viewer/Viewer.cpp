#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/surface_mesh.h"
#include "ExteriorPoissonRecon/ExteriorPoisson.h"

#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/Ply.h"
#include "Include/Hat.h"
#include "Misha/PlyVertexData.h"


const unsigned int Dim   = 3;
const unsigned int CoDim = 2;


// read from a text file
template< unsigned int Dim >
std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > ReadSamples( std::string fileName )
    {
	    std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > hermiteData;
	    if( Misha::GetFileExtension( fileName )==std::string( "ply" ) )
	{
		using Factory = VertexFactory::Factory< double, VertexFactory::PositionFactory< double, Dim >, VertexFactory::MatrixFactory< double , Dim , Dim > >;
		Factory factory = Factory( VertexFactory::PositionFactory< double , Dim >() , VertexFactory::MatrixFactory< double , Dim , Dim >( "skew" ) );
		std::vector< typename Factory::VertexType > vertices;
		std::vector< SimplexIndex< Dim-CoDim, unsigned int > > simplexIndices;
		int fileType;
		PLY::ReadSimplices( fileName , factory , vertices , simplexIndices , NULL , fileType );
		hermiteData.resize( vertices.size() );
		for( int i=0 ; i<vertices.size() ; i++ )
		{
			hermiteData[i].first = vertices[i].template get<0>();
			hermiteData[i].second = Hat::SkewSymmetricMatrix< double , Dim >( vertices[i].template get<1>() );
		}
	}
	else ERROR_OUT( "Only .ply files supported" );

	return hermiteData;
}

void oriented_circle(const Eigen::Vector3d &center, const Eigen::Vector3d &normal, double size, int n, Eigen::MatrixXd &C) {
    C.resize(n,3);
    Eigen::Vector3d s = Eigen::Vector3d::Random();
    s -= (s.dot(normal.normalized()) * normal.normalized());
    s.normalize();
    Eigen::Vector3d t = normal.cross(s).normalized();
    for (int i=0; i<n; i++) {
        float alpha =  ((float) i) / n * 2. * M_PI;
        C.row(i) = center + (s*cos(alpha) + t*sin(alpha)) * size;
    }
}

int main( int argc , char* argv[] )
{

    std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > slist;
    slist = ReadSamples<Dim>(argv[1]);

    std::vector<Eigen::VectorXd> tangents;
    // std::vector< Point<double,Dim>> points;
    Eigen::MatrixXd Points;

    int n=30;
    double s = 0.05;

    std::vector<Eigen::MatrixXd> disks;
    Points.resize(slist.size(),3);
    for (int i=0; i<slist.size(); i++){
            // points.push_back(slist[i].first);
            Points.row(i) = Eigen::RowVector3d(slist[i].first[0], slist[i].first[1], slist[i].first[2]);
            auto dual = slist[i].second;
            
            // conversion from wedge to cross product included 
            Eigen::Vector3d tangent = Eigen::Vector3d(dual[2], -dual[1], dual[0]);
            tangents.push_back(tangent);

            Eigen::MatrixXd C;
            oriented_circle(Eigen::Vector3d(slist[i].first[0], slist[i].first[1], slist[i].first[2]), tangent, s, n, C);
            disks.push_back(C);
    }

    // unite all disks in one polygon mesh
    Eigen::MatrixXd VD(disks.size()*n, 3);
    Eigen::MatrixXi FD(disks.size(), n);
    for (int d=0; d<disks.size(); d++) {
        VD.block(d*n,0,n,3) = disks[d]; 
        for (int t=0; t<n;t++){
            FD(d,t) = d*n+t;
        }
    }

    Eigen::MatrixXd face_colors(disks.size(), 3);
    for (int d=0; d<disks.size(); d++) {
        face_colors.row(d) = (tangents[d].normalized()/2. + Eigen::Vector3d(0.5, 0.5, 0.5));
    }
    
    // exterior product grid vis
    Eigen::RowVector3d corner = Points.colwise().minCoeff();
    double len = (Points.colwise().maxCoeff() - corner).maxCoeff() * 1.2;

    corner = Points.colwise().mean() - 0.5*len*Eigen::RowVector3d(1.,1.,1.); // / sqrt(2);

    int r = 15;
    double sigma = 0.1;
    double size = 1e-4;
    Eigen::MatrixXd GridCenters(r*r*r,3);
    std::vector<Eigen::MatrixXd> gdisks;
    for (int i=0; i<r; i++) {
        for (int j=0; j<r; j++) {
            for (int k=0; k<r; k++) {
                Eigen::RowVector3d c = Eigen::Vector3d(1./(r-1) * (int)i, 1./(r-1) * (int)j, 1./(r-1) * (int)k);
                c = corner + c*len;// *sqrt(2);
                GridCenters.row(i*r*r+j*r+k) = c;
                
                Hat::SkewSymmetricMatrix<double, Dim> s;
                for (int pi=0; pi<slist.size(); pi++) {
                    double weight = exp(- (c - Points.row(pi)).squaredNorm() / sigma);
                    s += weight * slist[pi].second;
                }
                Eigen::Vector3d tangent = Eigen::Vector3d(s[2], -s[1], s[0]);
                Eigen::MatrixXd C;
                std::cout << "tangent.norm: " << tangent.norm() << std::endl;
                oriented_circle(c.transpose(), tangent, tangent.norm()*size, n, C);
                gdisks.push_back(C);
            }
        }
    }

    // unite all grid disks
    Eigen::MatrixXd Vgd(gdisks.size()*n, 3);
    Eigen::MatrixXi Fgd(gdisks.size(), n);
    for (int d=0; d<gdisks.size(); d++) {
        Vgd.block(d*n,0,n,3) = gdisks[d]; 
        for (int t=0; t<n;t++){
            Fgd(d,t) = d*n+t;
        }
    }

    // -----------------------------
    // --------- POLYSCOPE ---------
    // -----------------------------

    polyscope::init();
    polyscope::view::upDir = polyscope::UpDir::ZUp;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::options::shadowBlurIters = 6;

    // auto grid = polyscope::registerPointCloud("Grid", GridCenters);
    auto gd = polyscope::registerSurfaceMesh("GridDisks", Vgd, Fgd);

    auto pc = polyscope::registerPointCloud("Input points", Points);
    pc->addVectorQuantity("Tangents", tangents);

    auto ds = polyscope::registerSurfaceMesh("Disks", VD, FD);
    ds->addFaceColorQuantity("Orientation", face_colors)->setEnabled(false);

    // update_visualization(P, T, EC, TF, level_vis[level], normalize_vectors, true);
    // polyscope::state::userCallback = myCallback;

    polyscope::show();

	return EXIT_SUCCESS;
}
