#pragma once

#include <iostream>
#include <stdexcept>

class Logger{
    public:
        Logger(std::ostream& os):os_(os){}
        template<typename T>
        Logger& operator<<(const T& t){
            os_ << t;
            return *this;
        }
        Logger& operator<<(std::ostream& (*pf)(std::ostream&)){
            os_ << pf;
            return *this;
        }
        Logger& printLog(const std::string& msg);
        Logger& printError(const std::string& msg);
        Logger& printWarning(const std::string& msg);
        void enable();
    private:
        std::ostream& os_;
        bool enabled_ = true;
};