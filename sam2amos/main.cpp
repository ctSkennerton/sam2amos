#include "api/BamMultiReader.h"
#include "foundation_AMOS.hh"
#include "amp.hh"
#include "StlExt.h"
#include "kseq.h"
#include <algorithm>
#include <string>
#include <map>
#include <utility>
#include <vector>
#include <cstring>
#include <getopt.h>

//using namespace BamTools;
//using namespace AMOS;
//using namespace std;
char comp_tab[] = {
    0,   1,	2,	 3,	  4,   5,	6,	 7,	  8,   9,  10,	11,	 12,  13,  14,	15,
    16,  17,  18,	19,	 20,  21,  22,	23,	 24,  25,  26,	27,	 28,  29,  30,	31,
    32,  33,  34,	35,	 36,  37,  38,	39,	 40,  41,  42,	43,	 44,  45,  46,	47,
    48,  49,  50,	51,	 52,  53,  54,	55,	 56,  57,  58,	59,	 60,  61,  62,	63,
    64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
	'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',	91,	 92,  93,  94,	95,
    64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
	'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};

std::string reverseComplement(std::string& str)
{
    
	int l = static_cast<int>(str.length());
    char * revcomp_str = new char[l+1];
    for (int i = 0; i <=l; i++) {
		revcomp_str[i]='\0';
    }
    int i, c0, c1;
    for (i = 0; i < l>>1; ++i) 
    {
        c0 = comp_tab[(int)str[i]];
        c1 = comp_tab[(int)str[l - 1 - i]];
        revcomp_str[i] = c1;
        revcomp_str[l - 1 - i] = c0;
    }
    if (l&1) 
    {
        revcomp_str[l>>1] = comp_tab[(int)str[l>>1]];
    }
    
    std::string ret = revcomp_str;
    delete [] revcomp_str;
    
    return ret;
}

std::string normalizeName(std::string& str) {
    // get the last two characters and see if they end in .f or .r
    std::string sub = str.substr(str.length() - 2);
    //std::cout <<sub<<std::endl;
    if (sub == ".f" || sub == ".r") {
        return str.substr(0, str.length() - 2);
    }
    return str;
}

void usage() {
    std::cout<<"sam2amos [options] <file.fa> <file.bam>...\n";
    std::cout<<"-h          This help message\n";
    std::cout<<"-b FILE     Output to bank file.  Default is to print\n";
    std::cout<<"            amos messages to screen\n";
    std::cout<<"-f INT      flags that must be set to output alignment\n";
    std::cout<<"            ALL bits must be set in alignment. default: 0\n";
    std::cout<<"-F INT      flags that must not be set to output alignment\n";
    std::cout<<"            if ANY bits are set alignment will not be outputted\n";
    std::cout<<"            Default: 0x900\n";
    std::cout<<"-a INT      Average insert size for the library when \n";
    std::cout<<"            the bam file contains paired reads\n";
    std::cout<<"-s INT      The standard deviation of the insert size\n";
    std::cout<<"            for the library"<<std::endl;
}

int processOptions(int argc, char *argv[], int& ave, int& stdev, std::string& bankFile, int& rejectFlag, int& acceptFlag)
{
    int c;
    int index;
    static struct option long_options [] = {
        {"help", no_argument, NULL, 'h'},
        {"bank",required_argument,NULL,'b'},
        {"average-insert",required_argument,NULL,'a'},
        {"stdev-insert",required_argument,NULL,'s'},
        {"reject-flag",required_argument,NULL,'F'},
        {"accept-flag",required_argument,NULL,'f'}
    };

    while( (c = getopt_long(argc, argv, "F:f:hb:a:s:", long_options, &index)) != -1 ) 
    {
        switch(c) 
        {
            case 'a':
                ave = atoi(optarg);
                break;
            case 'b':
                bankFile = optarg;
                break;
            case 'f':
                acceptFlag = (int)strtol(optarg, NULL, 0);
                break;
            case 'F':
                rejectFlag = (int)strtol(optarg, NULL, 0);
                break;
            case 's':
                stdev = atoi(optarg);
                break;
          case 'h':
          default:
                usage();
                exit(1);
        }
    }
    return optind;
}

int main(int argc, char ** argv) {

    
    AMOS::Message_t msg;
    
    // create stream objects for each of our message types
    AMOS::BankStream_t contig_bank(AMOS::Contig_t::NCODE);
    AMOS::BankStream_t read_bank(AMOS::Read_t::NCODE);
    AMOS::BankStream_t library_bank(AMOS::Library_t::NCODE);
    AMOS::BankStream_t fragment_bank(AMOS::Fragment_t::NCODE);
    
    if (argc < 3) {
        usage();
        return EXIT_FAILURE;
    }
    int library_mean = 0, library_stdev = 0, accept_flag = 0, reject_flag = 2304;
    std::string bank_name;
    int opt_idx = processOptions(argc, argv, library_mean, library_stdev, bank_name, reject_flag, accept_flag);
    if(opt_idx >= argc) {
        usage();
        return EXIT_FAILURE;
    }
    std::string fa_file = argv[opt_idx];
    if(opt_idx >= argc) {
        usage();
        return EXIT_FAILURE;
    }
    opt_idx++;

    std::vector<std::string> bam_files;
    while(opt_idx < argc) {
        bam_files.push_back(argv[opt_idx]);
        opt_idx++;
    }

    if(accept_flag & reject_flag) {
        std::cout<<"Bits are set for both accept and reject"<<std::endl;
        usage();
        return EXIT_FAILURE;
    }

    bool printmsg = bank_name.empty() ? true : false; 
    
    AMOS::Read_t read; 
    AMOS::Contig_t contig;
    AMOS::Tile_t tile;
    std::vector<AMOS::Tile_t> tiling;
    AMOS::Fragment_t fragment;
    
    try
    {
        if (!printmsg)
        {
            // Create the banks, alternatively use bank.open() to open an existing bank
            read_bank.create(bank_name);
            contig_bank.create(bank_name);
            library_bank.create(bank_name);
            fragment_bank.create(bank_name);
        }
        AMOS::ID_t read_id = 0;
        
        // provide some input & output filenames
        // attempt to open our BamMultiReader
        BamTools::BamMultiReader reader;
        if ( !reader.Open(bam_files) ) {
            std::cerr << "Could not open input BAM files." << std::endl;
            return 1; 
        }
        // retrieve 'metadata' from BAM files, these are required by BamWriter
        const BamTools::SamHeader header = reader.GetHeader();
                
        AMOS::Library_t library;
        AMOS::Distribution_t dist;

        // store the library iid and the @RG for later lookup
        typedef std::map<std::string, AMOS::ID_t> LibraryMap_t;  
        LibraryMap_t library_map;
        
        AMOS::ID_t library_id;
        // check to see if there are any read groups (@RG) in the header
        // this will correspond to the LIB messages as output
        if (header.HasReadGroups()) {
            BamTools::SamReadGroupConstIterator iter;
            
            for ( iter = header.ReadGroups.Begin(), library.clear(), library_id = 0; 
                 iter !=  header.ReadGroups.End(); 
                 iter++, library_id++) {
                
                library.setIID(library_id);
                library.setEID(iter->ID);
                
                library_map.insert(std::make_pair(iter->ID, library_id));
                
                if (iter->HasDescription()) {
                    library.setComment(iter->Description);
                }
                if (iter->HasPredictedInsertSize()) {
                    from_string(dist.mean, iter->PredictedInsertSize, std::dec);
                    library.setDistribution(dist);
                }
                if (printmsg) { library.writeMessage(msg); msg.write(std::cout); }
                else          { library_bank << library; }
            }
        } else {
            // we still need to make a libray
            library.setIID(1);
            if(library_mean) {
                dist.mean = library_mean;
            }
            if(library_stdev) {
                dist.sd = library_stdev;
            }
            library.setDistribution(dist);
            if (printmsg) { library.writeMessage(msg); msg.write(std::cout); }
            else          { library_bank << library; }
        }
        
        const BamTools::RefVector references = reader.GetReferenceData();

        
        // stores the name as string as the key and the internal
        // integer ID and the value for future lookup
        typedef std::map<std::string, int32_t> ContigMapID_t;
        ContigMapID_t contig_map;

        // read name base => fragment_t, ID_t
        typedef std::map<std::string, std::pair<AMOS::Fragment_t, AMOS::ID_t> > FragmentMap_t;
        FragmentMap_t fragment_map;
        
        // contig id => vector[tile_t]
        typedef std::map<int32_t, std::vector<AMOS::Tile_t> > TilingMap_t;
        TilingMap_t tiling_map;
        
        AMOS::ID_t fragment_id = 0;
        
        // iterate through all alignments, creating read, fragment and tiling
        // messages as appropriate
        BamTools::BamAlignment al;
        while ( reader.GetNextAlignmentCore(al) ) {
            if(al.AlignmentFlag & accept_flag != accept_flag) {
                continue;
            }
            
            if(al.AlignmentFlag & reject_flag) {
                continue;
            }

            al.BuildCharData();

            if (al.IsReverseStrand()) {
                al.QueryBases = reverseComplement(al.QueryBases);
                std::reverse(al.Qualities.begin(),al.Qualities.end());
            }
            // Internal IDs must be numeric (IID)
            read.setIID(++read_id);
            read.setSequence(al.QueryBases, al.Qualities);
            read.setType('E');
            
            // if the read is paired then we need to add
            // this read to a fragment, if one exists,
            // or make a new fragment
            if (al.IsPaired()) {
                std::string name;
                if(al.IsFirstMate()) {
                    name = al.Name + "/1";
                } else {
                    name = al.Name + "/2";
                }
                read.setEID(name);
                //std::string normalized_name = normalizeName(al.Name);
                FragmentMap_t::iterator fm_iter = fragment_map.find(al.Name);
                if (fm_iter != fragment_map.end()) {
                    // one of the fragments already exists
                    // get the name of the other one from the map
                    fm_iter->second.first.setReads(std::make_pair(read_id, fm_iter->second.second));
                    
                    // print out the fragment msg
                    if (printmsg) { fm_iter->second.first.writeMessage(msg); msg.write(std::cout); }
                    else          { fragment_bank << fm_iter->second.first; }
                    
                    // clean up since we-re done with this fragment
                    fragment_map.erase(fm_iter);
                } else {
                    // we've got to make a fragment
                    AMOS::Fragment_t frg;
                    frg.setIID(++fragment_id);
                    frg.setType('I');
                    frg.setEID(al.Name);
                    
                    // the libraries relate to the @RG tag
                    // no @RG, set to one
                    if (al.HasTag("RG")) {
                        std::string tag_data;
                        al.GetTag("RG", tag_data);
                        frg.setLibrary(library_map[tag_data]);
                    } else {
                        frg.setLibrary(1);
                    }
                    
                    fragment_map.insert(std::make_pair(al.Name, std::make_pair(frg,read_id)));
                }
                
            } else {
                read.setEID(al.Name);
            }
            // print out the read msg
            if (printmsg) { read.writeMessage(msg); msg.write(std::cout); }
            else          { read_bank << read; }
            
            // make the tilling for the read
            tile.source = read_id;
            tile.offset = al.Position;
            AMOS::Pos_t qstart = al.Position, qend = al.Position + al.Length;
            
            // take into account the soft cliping
            if (al.CigarData.front().Type == 'S') {
                qstart += al.CigarData.front().Length; 
            }
            if (al.CigarData.back().Type == 'S') {
                qend -= al.CigarData.back().Length;
            }
            
            // account for reverse complements
            AMOS::Range_t tile_range;
            if (al.IsReverseStrand()) {
                tile_range.begin = qend;
                tile_range.end = qstart;
            
            } else {
                tile_range.begin = qstart;
                tile_range.end = qend;
            }
            tile.range = tile_range;
            
            // add the tile to the contig
            tiling_map[al.RefID].push_back(tile);
        
        }
        
        // now go through the fasta file and get the contig sequences
        gzFile fp = gzopen(fa_file.c_str(), "rb");
        kseq_t * seq;
        
        seq = kseq_init(fp);
        while (kseq_read(seq) >= 0 ) 
        {
            ContigMapID_t::iterator ctg_iter;
            ctg_iter = contig_map.find(seq->name.s);
            
            if(ctg_iter == contig_map.end()) {
                contig.clear();
                contig.setIID(ctg_iter->second);
                contig.setEID(seq->name.s);
                contig.setSequence(seq->seq.s, seq->seq.s);
                contig.setReadTiling(tiling_map[ctg_iter->second]);
                // print out the contig msg
                if (printmsg) { contig.writeMessage(msg); msg.write(std::cout); }
                else          { contig_bank << contig; }
                
                // clean up since we're done with this contig
                contig_map.erase(ctg_iter);
            }
        }
                
        if (!printmsg)
        {
            // You should now be able to open the bank with bankViewer
            // and see the multiple alignment described above.
            read_bank.close();
            contig_bank.close();
        }
    }
    catch (AMOS::Exception_t & e)
    {
        std::cerr << "ERROR: -- Fatal AMOS Exception --\n" << e;
        return EXIT_FAILURE;
    }
    catch (std::exception& e) {
        std::cerr << e.what()<<std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
