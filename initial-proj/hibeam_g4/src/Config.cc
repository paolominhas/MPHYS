// -------------------------------------------------------
// Implementation of the Config class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------------

#include "Config.hh"

#include "G4SystemOfUnits.hh"

#include <fstream>
#include <getopt.h>
#include <sstream>

static struct option optlong[] = {
	{"help", 0, NULL, 'h'}, {"gui", 0, NULL, 'g'},   {"mac", 1, NULL, 'm'},
	{"seed", 1, NULL, 's'}, {"input", 1, NULL, 'i'}, {"config", 1, NULL, 'c'}, {"threads", 1, NULL, 't'},
};
Config::Config() { SetDefault(); }

Config::Config(int argc, char** argv)
{
	SetDefault();
	status = ParseCmd(argc, argv);
}

Config::~Config() {}

void Config::ParseLine(std::stringstream& lineStream, std::vector<std::string>& res)
{
	std::string lastStr;
	lineStream >> lastStr;
	if(lastStr.empty())
		return;
	else
	{
		res.emplace_back(lastStr);
		return ParseLine(lineStream, res);
	}
}

void Config::ParseConfig(const std::string& namefile)
{
	std::cout << "start reading" << std::endl;
	std::ifstream ifs(namefile.c_str());
	if(ifs.is_open())
	{
		const std::string CommentSymbol("#");

		std::string temp_line;
		while(std::getline(ifs, temp_line))
		{
			std::stringstream stream(temp_line);
			std::string testComment(stream.str());
			auto it_comment = testComment.find(CommentSymbol);
			if(it_comment != std::string::npos)
			{
				// std::cout<<"!> Skip "<<test<<std::endl;
				continue;
			}

			std::vector<std::string> out_list;
			ParseLine(stream, out_list);

			switch(out_list.size())
			{
				case 2:
					par[out_list[0]] = out_list[1];
					break;
				default:
					std::cout << "!> stream of line from config file with unknown parsing " << out_list.size() << " \n";
					break;
			}

		}
	}
}

void Config::ParsePreviousParams(std::map<std::string,double>* params)
{
	
	auto temp_pair = params->find("MCPL_Inputfile");
	if(temp_pair != params->end())
	{
		par["MCPL_Inputfile"] = temp_pair->second;
	}
	
}

int Config::ParseCmd(int argc, char** argv)
{

	auto print_help = [&argv]() {
		std::cout << "!> Wrong number of parameters!\n";
		std::cout << "!> Example of use:\n";
		std::cout << "!> " << argv[0];
		std::cout << "!> [-h] [--help] [-g] [--gui] [-s seed] [--seed seed] [-i inputfile] [--input inputfile] [-m "
			"run.mac] [--mac run.mac] [-t threads] [--threads threads] [-c config.par] [--config config.par] Outputfile.root \n";
		std::cout << " \n";
	};

	if(argc < 2)
	{
		print_help();
		return -1;
	}
	int option_char;
	std::string nameI, nameM, nameC("testconfig.par");
	std::string seed, threads;
	while((option_char = getopt_long(argc, argv, "+hgst:m:i:c:", optlong, NULL)) != EOF)
		switch(option_char)
		{
			case 'h':
				print_help();
				return -1;
				break;
			case 'g':
				std::cout << "Gui mode " << std::endl;
				par["Gui"] = "1";
				break;
			case 's':
				std::cout << "Seed for Parallel runs :" << optarg << std::endl;
				seed = optarg;
				par["HEPRand_Seed"] = seed;
				break;
			case 'i':
				std::cout << "Inputfile of Event :" << optarg << std::endl;
				nameI = optarg;
				par["InputFile"] = nameI;
				break;
			case 'm':
				std::cout << "Macro File :" << optarg << std::endl;
				nameM = optarg;
				par["MacroFile"] = nameM;
				break;
				case 'c':
				std::cout << "Configuration File :" << optarg << std::endl;
				nameC = optarg;
				par["Config"] = nameC;
				break;
			case 't':
				std::cout << "Number of threads :" << optarg << std::endl;
				threads = optarg;
				par["Threads"] = threads;
				break;
			case '?':
				print_help();
				return -1;
		}

	std::string name_out;

	if(optind == argc)
	{
		print_help();
		return -1;
	}
	else
	{
		name_out = argv[optind];
		par["Output_Namefile"] = name_out;
	}

	ParseConfig(nameC);
	
	return 0;
}

void Config::display()
{
	for (const auto& [key, value] : par)
        std::cout << key << " -> " << value << std::endl;

}

void Config::CheckConfig(){ 
	display(); 
}

void Config::SetDefault()
{
	par["Gui"] = "0";
	par["Threads"] = "1";
	par["Output_Namefile"] = "Default_Output.root";
	par["Geometry_Namefile"] = "geometry.root";
	par["CheckOverlaps"] = "0";

	par["Source"] = "gps";
	par["MCPL_Inputfile"] = "none";
	par["Detectors"] = "";
	par["Sampling_Detectors"] = "";
	par["WriteHistograms"] = "1";
	par["WriteTree"] = "0";
}

bool Config::IsAvailable(const std::string& key)
{
	return (par.find(key) != par.end());
}

std::string Config::GetString(const std::string& key)
{
	std::string ret = "";
	if(IsAvailable(key)) ret = par[key];
	else std::cout << "Parameter " << key << " unavailable!" << std::endl;
	return ret;
}

int Config::GetInt(const std::string& key)
{
	int ret = 0;
	if(IsAvailable(key)) ret = std::stoi(par[key]);
	else std::cout << "Parameter " << key << " unavailable!" << std::endl;
	return ret;
}

double Config::GetDouble(const std::string& key)
{
	double ret = 0.;
	if(IsAvailable(key)) ret = std::stod(par[key]);
	else std::cout << "Parameter " << key << " unavailable!" << std::endl;
	return ret;
}
