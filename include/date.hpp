#pragma once
#include<chrono>
#include<string>
#include<sstream>

inline const std::string getCurrentSystemTime();
// 获取年月日,用_隔开
inline const std::string getDate();

const std::string getCurrentSystemTime(){
	auto tt = std::chrono::system_clock::to_time_t
	(std::chrono::system_clock::now());
	struct tm* ptm = localtime(&tt);
	char date[60] = {0};
	sprintf(date, "%d-%02d-%02d %02d:%02d:%02d",
		(int)ptm->tm_year + 1900,(int)ptm->tm_mon + 1,(int)ptm->tm_mday,
		(int)ptm->tm_hour,(int)ptm->tm_min,(int)ptm->tm_sec);
	return std::string(date);
}

// 获取年月日,用_隔开
const std::string getDate(){
	return getCurrentSystemTime().substr(0,4) + "_" + getCurrentSystemTime().substr(5,2) + "_" + getCurrentSystemTime().substr(8,2);
}