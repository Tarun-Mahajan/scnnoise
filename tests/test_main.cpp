// #define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include <catch2/catch_all.hpp>
#include <catch2/internal/catch_clara.hpp>
#include <string>
#include "../src/scnnoise/utils.hpp"
#include "test_GRN_global_variables.hpp"

std::string GRN_filepath = "temp";
unsigned int num_nodes = 10;

int main( int argc, char* argv[] )
{
  Catch::Session session; // There must be exactly one instance

  // int height = 0; // Some user variable you want to be able to set
  // Build a new parser on top of Catch's
  using namespace Catch::Clara;
  auto cli
    = session.cli() // Get Catch's composite command line parser
    | Opt( GRN_filepath, "GRN_filepath" ) // bind variable to a new option, with a hint string
        ["-g"]["--grn_file"]    // the option names it will respond to
        ("Path to GRN file.")
    | Opt( num_nodes, "num_nodes" ) // bind variable to a new option, with a hint string
        ["-n"]["--num_nodes"]    // the option names it will respond to
        ("Number of genes in the GRN.");        // description string for the help output

  // Now pass the new composite back to Catch so it uses that
  session.cli( cli );

  // Let Catch (using Clara) parse the command line
  int returnCode = session.applyCommandLine( argc, argv );
  if( returnCode != 0 ) // Indicates a command line error
      return returnCode;

  // if set on the command line then 'height' is now set at this point
  // if( height > 0 )
  //     std::cout << "height: " << height << std::endl;

  return session.run();
}
