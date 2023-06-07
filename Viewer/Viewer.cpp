#include "polyscope/polyscope.h"
#include "ExteriorPoissonRecon/ExteriorPoisson.h"

#include "Misha/Miscellany.h"
#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/Ply.h"
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

int main( int argc , char* argv[] )
{


    std::vector< std::pair< Point< double , Dim > , Hat::SkewSymmetricMatrix< double , Dim > > > slist;
    slist = ReadSamples<Dim>(argv[1]);

    // -----------------------------
    // --------- POLYSCOPE ---------
    // -----------------------------

    polyscope::init();
    polyscope::view::upDir = polyscope::UpDir::ZUp;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::options::shadowBlurIters = 6;

    // update_visualization(P, T, EC, TF, level_vis[level], normalize_vectors, true);
    // polyscope::state::userCallback = myCallback;

    polyscope::show();

	return EXIT_SUCCESS;
}
