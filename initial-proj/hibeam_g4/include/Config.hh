// -------------------------------------------------
// Definition of the Config class
// Created by C.Rappold (c.rappold@gsi.de)
//--------------------------------------------------

#ifndef CONFIG_h
#define CONFIG_h

#include <iostream>
#include <string>
#include <vector>
#include <map>

class Config
{
	public:
		Config();
		Config(int argc, char** argv);
		~Config();
		void ParseConfig(const std::string& namefile);
		void ParsePreviousParams(std::map<std::string,double>* params);
		int ProperConf() { return status; }

		bool IsAvailable(const std::string& key);
		std::string GetString(const std::string& key);
		int GetInt(const std::string& key);
		double GetDouble(const std::string& key);
		
		void CheckConfig();

	private:
		std::map<std::string,std::string> par;
		int status;

		int ParseCmd(int argc, char** argv);
		void ParseLine(std::stringstream& lineStream, std::vector<std::string>& res);
		void display();
		void SetDefault();
};

#endif
