#include "Logger.hxx"


Logger& Logger::printLog(const std::string& msg){
    if(enabled_){
        os_ << "[LOG] "<< msg;
    }
    return *this;
}

Logger& Logger::printError(const std::string& msg){
    if(enabled_){
        os_ << "[ERROR] "<< msg;
    }
    return *this;
}

Logger& Logger::printWarning(const std::string& msg){
    if(enabled_){
        os_ << "[WARNING] "<< msg;
    }
    return *this;
}

void Logger::enable(){
    enabled_ = true;
}