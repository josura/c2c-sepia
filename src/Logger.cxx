#include "Logger.hxx"


Logger& Logger::printLog(const std::string& msg){
    if(enabled_){
        os_ << "[LOG] "<< msg;
    }
    return *this;
}

void Logger::enable(){
    enabled_ = true;
}