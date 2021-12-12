// vim: set sw=4 ts=4 sts=4 expandtab :
/* -------------------------------------------------------------------------- *
 *                           MMB (MacroMoleculeBuilder)                       *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * Copyright (c) 2011-12 by the Author.                                       *
 * Author: Samuel Flores                                                      *
 *                                                                            *
 * See RNABuilder.cpp for the copyright and usage agreement.                  *
 * -------------------------------------------------------------------------- */

#include "Utils.h"
#include <MMBLogger.h>

#include <cstring>
#include <functional>
#include <iostream>
#include <sstream>

#ifdef LEPTON_ENABLED
#include <Lepton.h>
#endif

#ifdef _WINDOWS
#include <windows.h>
#include <memory>
#else
    #include <fcntl.h>
    #include <sys/utsname.h>
    #ifdef HAVE_COPY_FILE_RANGE
        #ifndef _GNU_SOURCE
            #define _GNU_SOURCE
        #endif // _GNU_SOURCE
        #include <unistd.h>
    #endif // HAVE_COPY_FILE_RANGE
    #ifdef HAVE_SENDFILE
        #include <sys/sendfile.h>
    #endif // HAVE_SENDFILE
#endif // _WINDOWS

using namespace std;
using namespace SimTK;

#ifdef _WINDOWS

class SecurityDescriptorWrapper {
public:
    SecurityDescriptorWrapper(DWORD bufSize) {
        if (bufSize == 0)
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Invalid security descriptor buffer size");

        m_buffer = std::make_unique<char *>(new char[bufSize]);
    }

    PSECURITY_DESCRIPTOR get() {
        return static_cast<PSECURITY_DESCRIPTOR>(m_buffer.get());
    }

private:
    std::unique_ptr<char *> m_buffer;
};

inline
std::string getErrorString(DWORD error) {
	LPSTR buf = nullptr;

    if (FormatMessageA(
            FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
		    NULL,
            error,
            MAKELANGID(LANG_NEUTRAL, SUBLANG_NEUTRAL),
            (LPSTR)&buf,
            0,
            nullptr
	    ) == 0) {
        return "Unknown error";
    }

    std::string msg(buf);
    LocalFree(buf);

    return msg;
}

inline
bool hasAccessRight(LPCSTR path, DWORD genericAccessRights) {
    SECURITY_DESCRIPTOR secDesc;

    DWORD len = 0;
    if (GetFileSecurityA(path, OWNER_SECURITY_INFORMATION | GROUP_SECURITY_INFORMATION | DACL_SECURITY_INFORMATION, NULL, NULL, &len) == FALSE) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Failed to get security information for file " << path << std::endl);
    }

    SecurityDescriptorWrapper sd{ len };
    if (GetFileSecurityA(path, OWNER_SECURITY_INFORMATION | GROUP_SECURITY_INFORMATION | DACL_SECURITY_INFORMATION, sd.get(), len, &len) == FALSE) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Failed to get security information for file " << path << std::endl);
    }

    HANDLE hToken = NULL;
    HANDLE hImpersonatedToken = NULL;
    if (OpenProcessToken(GetCurrentProcess(), TOKEN_IMPERSONATE | TOKEN_QUERY | TOKEN_DUPLICATE | STANDARD_RIGHTS_READ, &hToken) == FALSE) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Failed to get security information for file " << path << std::endl);
    }

    if (DuplicateToken(hToken, SecurityImpersonation, &hImpersonatedToken) == FALSE) {
        CloseHandle(hToken);
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Failed to get security information for file " << path << std::endl);
    }

    GENERIC_MAPPING mapping = { 0xFFFFFFFF };
    PRIVILEGE_SET privileges = { 0 };
    DWORD grantedAccess = 0;
    DWORD privilegesLen = sizeof(privileges);

    mapping.GenericRead = FILE_GENERIC_READ;
    mapping.GenericWrite = FILE_GENERIC_WRITE;
    mapping.GenericExecute = FILE_GENERIC_EXECUTE;
    mapping.GenericAll = FILE_ALL_ACCESS;
    BOOL result = FALSE;

    MapGenericMask(&genericAccessRights, &mapping);
    auto success = AccessCheck(sd.get(), hImpersonatedToken, genericAccessRights, &mapping, &privileges, &privilegesLen, &grantedAccess, &result);
    CloseHandle(hImpersonatedToken);
    CloseHandle(hToken);

    if (success == FALSE)
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Failed to get security information for file " << path << std::endl);

    return result == TRUE;
}

inline
bool hasReadAccess(LPCSTR path) {
    return hasAccessRight(path, GENERIC_READ);
}

inline
bool hasWriteAccess(LPCSTR path) {
    return hasAccessRight(path, GENERIC_READ | GENERIC_WRITE);
}

inline
std::wstring multibyteToWideString(const std::string &str) {
    auto ret = MultiByteToWideChar(
        CP_UTF8,
        MB_ERR_INVALID_CHARS,
        str.c_str(),
        -1,
        NULL,
        0
    );

    if (ret < 1) {
        throw new std::runtime_error("Cannot calculate size of the array for widechar string");
    }

    auto buf = std::make_unique<WCHAR[]>(ret);

    ret = MultiByteToWideChar(
        CP_UTF8,
        MB_ERR_INVALID_CHARS,
        str.c_str(),
        -1,
        buf.get(),
        ret
    );

    if (ret < 1) {
        throw new std::runtime_error("Cannot convert to widechar string");
    }

    return std::wstring(buf.get());
}

inline
std::wstring unicodeFullPath(const std::string &path) {
    auto wpath = multibyteToWideString(path);

    auto buf = std::make_unique<WCHAR[]>(32768);
    if (GetFullPathNameW(wpath.c_str(), 32768, buf.get(), nullptr) == 0) {
        throw new std::runtime_error("Cannot get full path name");
    }

    return std::wstring(buf.get());
}

#endif // _WINDOWS

int myMkdir(const std::string & directoryPath) {
#ifdef _WINDOWS
    MMBLOG_FILE_FUNC_LINE(INFO, " You are asking to create the directory  " << directoryPath << std::endl);

    WIN32_FIND_DATA fData;
    HANDLE hDir = FindFirstFileA(directoryPath.c_str(), &fData);
    if (hDir == INVALID_HANDLE_VALUE) {
        if (CreateDirectoryA(directoryPath.c_str(), NULL) == FALSE) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, " Failed to create directory " << directoryPath << ", error " << GetLastError() << std::endl);
        }
    }
    hDir = FindFirstFileA(directoryPath.c_str(), &fData);

    if (!hasWriteAccess(directoryPath.c_str())) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Heads up! Found that we do NOT have write access to a directory called " << directoryPath << "  " << std::endl);
        CloseHandle(hDir);
        return 1;
    }

    return 0;
#else
    MMBLOG_FILE_FUNC_LINE(INFO, " You are asking to create the directory  "<<directoryPath<<std::endl);
    if (!(opendir(directoryPath.c_str()))){
        MMBLOG_FILE_FUNC_LINE(INFO, " opendir failed to open directory "<<directoryPath<<" . Will now create this directory.  " <<std::endl);
        const int dir_err = mkdir(directoryPath.c_str(), S_IRWXU | S_IRGRP | S_IXGRP | S_IROTH | S_IXOTH);
        MMBLOG_FILE_FUNC_LINE(INFO, ": mkdir returned : "<<dir_err <<std::endl);
        if (-1 == dir_err)
        {
            printf("Error creating directory!n");
            MMBLOG_FILE_FUNC_LINE(INFO, " Failed to  create directory "<<directoryPath<<"  " <<std::endl);
            exit(1);
        } else if (0 == dir_err) {
            MMBLOG_FILE_FUNC_LINE(INFO, " Successfully created directory "<<directoryPath<<" " <<std::endl);
        } else {
            MMBLOG_FILE_FUNC_LINE(INFO, " An unexpected error occurred when creating the directory "<<directoryPath<<" " <<std::endl);
        }
    }
    if (access((directoryPath ).c_str(), R_OK) == 0) {
        MMBLOG_FILE_FUNC_LINE(INFO, " Found that we have read access to a directory called "<< directoryPath   <<" . So far so good."<<std::endl);
        return 0;
    } else {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Heads up! Found that we do NOT have read access to a directory called "<<  directoryPath  <<"  "<<std::endl);
        return 1;
    }
#endif // _WINDOWS
}

int myChdir(const std::string & directoryPath) {
    MMBLOG_FILE_FUNC_LINE(INFO, " About to attempt changing directory to " << directoryPath << " . " << std::endl);
#ifdef _WINDOWS
    if (SetCurrentDirectoryA(directoryPath.c_str()) == FALSE) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Unable to change directory to " << directoryPath << " . Exiting now." << std::endl);
    }
    MMBLOG_FILE_FUNC_LINE(INFO, " Was able to successfully change directory to " << directoryPath << " . " << std::endl);
    return 0;
#else
    if (chdir(directoryPath.c_str())){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Unable to change directory to "<<directoryPath<<" . Exiting now."<<std::endl);
        return 1;
    } else {
        MMBLOG_FILE_FUNC_LINE(INFO, " Was able to successfully change directory to "<<directoryPath<<" . "<<std::endl);
        return 0;
    }
#endif // _WINDOWS
}

/* Copy file function variants */
#ifdef _WINDOWS

static const std::wstring LONG_PATH_LIMIT_PREFIX = L"\\\\?\\";
CopyFileResult mmbCopyFile(const std::string &sourceFileName, const std::string &destinationFileName) {
    try {
        auto src = LONG_PATH_LIMIT_PREFIX + unicodeFullPath(sourceFileName);
        auto dst = LONG_PATH_LIMIT_PREFIX + unicodeFullPath(destinationFileName);

        if (!CopyFileW(src.c_str(), dst.c_str(), FALSE)) {
            auto err = GetLastError();
            if (err == ERROR_DISK_FULL) {
                return CopyFileResult::No_Space;
            }
            MMBLOG_PLAIN(WARNING, "IO error while copying a file: " << getErrorString(err) << "\n");
            return CopyFileResult::Io_Error;
        }
        return CopyFileResult::Success;
    } catch (const std::runtime_error &ex) {
        MMBLOG_PLAIN(WARNING, ex.what() << "\n");
        return CopyFileResult::Io_Error;
    }
}

#else

static
auto getFileSizeFromFd(int fd) {
    struct stat st;

    int ret = fstat(fd, &st);
    if (ret == -1) {
        return -1L;
    }

    return st.st_size;
}

#ifdef HAVE_COPY_FILE_RANGE

static
CopyFileResult mmbCopyFile_copy_file_range(const std::string &sourceFileName, const std::string &destinationFileName) {
    int in = open(sourceFileName.c_str(), O_RDONLY);
    if (in == -1) {
        MMBLOG_PLAIN(WARNING, "Cannot open file " << sourceFileName << " for reading");

        return CopyFileResult::Io_Error;
    }

    const auto srcSize = getFileSizeFromFd(in);
    if (srcSize == -1) {
        MMBLOG_PLAIN(WARNING, "Cannot determine size of the file to copy:" << strerror(errno));
        close(in);

        return CopyFileResult::Io_Error;
    }

    int out = creat(destinationFileName.c_str(), S_IWUSR | S_IRUSR);
    if (out == -1) {
        MMBLOG_PLAIN(WARNING, "Cannot open file " << destinationFileName << " for writing");
	    close(in);

        return CopyFileResult::Io_Error;
    }

    CopyFileResult retCode = CopyFileResult::Success;
    auto bytesToCopy = srcSize;
    while (bytesToCopy > 0) {
        auto ret = copy_file_range(in, nullptr, out, nullptr, bytesToCopy, 0);
        if (ret == -1) {
            if (errno == ENOSPC || errno == EDQUOT) {
                retCode = CopyFileResult::No_Space;
            } else {
                MMBLOG_PLAIN(WARNING, "IO error while copying a file: " << strerror(errno));

                retCode = CopyFileResult::Io_Error;
            }
            goto out;
        }

        bytesToCopy -= ret;
    }

out:
    close(in);
    close(out);

    return retCode;
}

#endif // HAVE_COPY_FILE_RANGE

CopyFileResult mmbCopyFile_direct(const std::string &sourceFileName, const std::string &destinationFileName) {
    const size_t BUFSIZE = 8 * 1024;
    auto buf = std::make_unique<char[]>(BUFSIZE);

    int in = open(sourceFileName.c_str(), O_RDONLY);
    if (in == -1) {
        MMBLOG_PLAIN(WARNING, "Cannot open file " << sourceFileName << " for reading");

        return CopyFileResult::Io_Error;
    }

    const auto srcSize = getFileSizeFromFd(in);
    if (srcSize == -1) {
        MMBLOG_PLAIN(WARNING, "Cannot determine size of the file to copy:" << strerror(errno));
        close(in);

        return CopyFileResult::Io_Error;
    }

    int out = creat(destinationFileName.c_str(), S_IWUSR | S_IRUSR);
    if (out == -1) {
        MMBLOG_PLAIN(WARNING, "Cannot open file " << destinationFileName << " for writing");
	    close(in);

        return CopyFileResult::Io_Error;
    }

    CopyFileResult retCode = CopyFileResult::Success;
    auto bytesToCopy = srcSize;
    while (bytesToCopy > 0) {
        auto blkBytesToCopy = read(in, buf.get(), BUFSIZE);
        if (blkBytesToCopy == -1) {
            MMBLOG_PLAIN(WARNING, "IO read error while copying a file: " << strerror(errno));

            goto out;
        }

        while (blkBytesToCopy > 0) {
            auto ret = write(out, buf.get(), blkBytesToCopy);
            if (ret == -1) {
                if (errno == ENOSPC || errno == EDQUOT) {
                    retCode = CopyFileResult::No_Space;
                } else {
                    MMBLOG_PLAIN(WARNING, "IO write error while copying a file: " << strerror(errno));

                    retCode = CopyFileResult::Io_Error;
                }

                goto out;
            }

            blkBytesToCopy -= ret;
            bytesToCopy -= ret;
        }
    }

    if (fsync(out) == -1) {
        if (errno == ENOSPC || errno == EDQUOT) {
            retCode = CopyFileResult::No_Space;
        } else {
            MMBLOG_PLAIN(WARNING, "IO write error while copying a file: " << strerror(errno));

            retCode = CopyFileResult::Io_Error;
        }
    }

out:
    close(in);
    close(out);

    return retCode;
}

#ifdef HAVE_SENDFILE

static
CopyFileResult mmbCopyFile_sendfile(const std::string &sourceFileName, const std::string &destinationFileName) {
    int in = open(sourceFileName.c_str(), O_RDONLY);
    if (in == -1) {
        MMBLOG_PLAIN(WARNING, "Cannot open file " << sourceFileName << " for reading");

        return CopyFileResult::Io_Error;
    }

    const auto srcSize = getFileSizeFromFd(in);
    if (srcSize == -1) {
        MMBLOG_PLAIN(WARNING, "Cannot determine size of the file to copy: " << strerror(errno));
        close(in);

        return CopyFileResult::Io_Error;
    }

    int out = creat(destinationFileName.c_str(), S_IWUSR | S_IRUSR);
    if (out == -1) {
        MMBLOG_PLAIN(WARNING, "Cannot open file " << destinationFileName << " for writing");
	    close(in);

        return CopyFileResult::Io_Error;
    }

    CopyFileResult retCode = CopyFileResult::Success;
    auto bytesToCopy = srcSize;
    while (bytesToCopy > 0) {
        auto ret = sendfile(out, in, nullptr, srcSize);
        if (ret == -1) {
            if (errno == ENOSPC || errno == EDQUOT) {
                retCode = CopyFileResult::No_Space;
            } else if ((errno == EINVAL || errno == ENOSYS) && bytesToCopy == srcSize) {
                /* Fall back to direct copying - system we are running on may not implement sendfile() the way we need.
                 * This makes sense only if the first call to sendfile() failes, subsequent failures with the error codes
                 * above indicate a problem somewhere else */
                close(in);
                close(out);
                return mmbCopyFile_direct(sourceFileName, destinationFileName);
            } else {
                MMBLOG_PLAIN(WARNING, "IO error while copying a file: " << strerror(errno));

                retCode = CopyFileResult::Io_Error;
            }
            goto out;
        }

        bytesToCopy -= ret;
    }

    if (fsync(out) == -1) {
        if (errno == ENOSPC || errno == EDQUOT) {
            retCode = CopyFileResult::No_Space;
        } else {
            MMBLOG_PLAIN(WARNING, "IO error while copying a file: " << strerror(errno));

            retCode = CopyFileResult::Io_Error;
        }
    }

out:
    close(in);
    close(out);

    return retCode;
}

static
std::tuple<bool, int, int, int> getLinuxKernelVersion() {
    struct utsname kernelInfo;

    if (uname(&kernelInfo) != 0) {
        return { false, 0, 0, 0 };
    }

    if (std::strcmp(kernelInfo.sysname, "Linux") != 0) {
        return { false, 0, 0, 0 };
    }

    std::vector<std::string> parts{};

    const size_t len = strlen(kernelInfo.release);
    size_t from = 0;

    for (size_t to = 1; to < len; to++) {
        if (kernelInfo.release[to] == '.') {
            parts.push_back(std::string(&kernelInfo.release[from], to - from));
            from = to + 1;
        }
    }

    if (parts.size() < 2)
        return { false, 0, 0, 0 };

    try {
        auto major = std::stoi(parts[0]);
        auto minor = std::stoi(parts[1]);

        if (major > 2) {
            return { true, major, minor, 0 };
        }

        // Check for 2.6.x naming scheme
        if (parts.size() < 3) {
            return { false, 0, 0, 0 };
        }

        auto minorTwo = std::stoi(parts[2]);
        return { true, major, minor, minorTwo };
    } catch (const std::invalid_argument &) {
        return { false, 0, 0, 0 };
    } catch (const std::out_of_range &) {
        return { false, 0, 0, 0 };
    }
}

static
std::function<CopyFileResult (const std::string &, const std::string &)> getMmbCopyFileImpl() {
    auto lnxKernVer = getLinuxKernelVersion();

    if (!std::get<0>(lnxKernVer)) {
        // We are not running on Linux or cannot determine the kernel version
        return mmbCopyFile_direct;
    }

    auto maj = std::get<1>(lnxKernVer);
    auto min = std::get<2>(lnxKernVer);
    auto minTwo = std::get<3>(lnxKernVer);

#ifdef HAVE_COPY_FILE_RANGE
    if (maj > 5) {
        return mmbCopyFile_copy_file_range;
    } else if (maj == 5 && min >= 3) {
        return mmbCopyFile_copy_file_range;
    }
#endif // HAVE_COPY_FILE_RANGE

#ifdef HAVE_SENDFILE
    if (maj >= 3 || (maj == 2 && min == 6 && minTwo >= 33)) {
        return mmbCopyFile_sendfile;
    }
#endif // HAVE_SENDFILE

    return mmbCopyFile_direct;
}

#endif // HAVE_SENDFILE

CopyFileResult mmbCopyFile(const std::string &sourceFileName, const std::string &destinationFileName) {
    static const std::function<CopyFileResult (const std::string &, const std::string &)> impl = getMmbCopyFileImpl();

    return impl(sourceFileName, destinationFileName);
}

#endif // _WINDOWS

void closingMessage() {
    std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<" Questions or problems? Please use the public forum at https://simtk.org/forums/viewforum.php?f=359&amp;sid=770b2f0d333ecb740d8c2f9e7e80e51c  .. or email me."<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<" Please support continued development. "<<std::endl;
    //std::cout<<__FILE__<<":"<<__LINE__<<" Suggested donation: "<<std::endl;
    //std::cout<<__FILE__<<":"<<__LINE__<<" Industry (per user): 1000 EUR "<<std::endl;
    //std::cout<<__FILE__<<":"<<__LINE__<<" By bank transfer to IBAN: SE0750000000053680279418 , SWIFT: ESSESESS "<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<" Industrial and other inquiries: samuel.flores@scilifelab.se "<<std::endl;
    std::cout<<__FILE__<<":"<<__LINE__<<std::endl;
};

CheckFile::CheckFile(const String & myFileName){
    fileName = myFileName;
#ifndef _WINDOWS
    stat(  fileName.c_str(), &st);
    MMBLOG_FILE_FUNC_LINE(INFO, " stat .st_mode= "<<st.st_mode <<std::endl);
#endif // _WINDOWS
}
bool CheckFile::isDirectory(){
    MMBLOG_FILE_FUNC_LINE(INFO, " stat S_IFDIR = "<<S_IFDIR <<std::endl);
    MMBLOG_FILE_FUNC_LINE(INFO, " stat .st_mode= "<<st.st_mode <<std::endl);
    MMBLOG_FILE_FUNC_LINE(INFO, " st.st_mode & S_IFMT "<< (st.st_mode & S_IFMT)  <<std::endl);
    //return (st.st_mode == S_IFDIR);	
    return ((st.st_mode & S_IFMT) == S_IFDIR);
}
bool CheckFile::ownerCanRead() {
#ifdef _WINDOWS
    return hasReadAccess(fileName.c_str());
#else
    MMBLOG_FILE_FUNC_LINE(INFO, " st.st_mode & S_IRUSR"<< (st.st_mode & S_IRUSR) <<std::endl);
    //return (st.st_mode == S_IFDIR);	
    return ((st.st_mode ) & S_IRUSR);
#endif // _WINDOWS
}
bool CheckFile::ownerCanWrite() {
#ifdef _WINDOWS
    return hasWriteAccess(fileName.c_str());
#else
    MMBLOG_FILE_FUNC_LINE(INFO, " st.st_mode & S_IWUSR"<< (st.st_mode & S_IWUSR) <<std::endl);
    //return (st.st_mode == S_IFDIR);	
    return ((st.st_mode ) & S_IWUSR);
#endif // _WINDOWS
}

void CheckFile::validateNonZeroSize(){
    MMBLOG_FILE_FUNC_LINE(INFO, " About to check that "<<fileName<<" has nonzero size.."<<std::endl);
    if ( st.st_size == 0){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Apparently "<<fileName<<" has size "<<st.st_size <<" . Dying now."<<std::endl);
    } else {
        MMBLOG_FILE_FUNC_LINE(INFO, " Apparently "<<fileName<<" has size "<<st.st_size <<" . This seems OK."<<std::endl);
    }
}

void CheckFile::validateExists(){
    if(stat(fileName.c_str(), &st) != 0){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " ERROR: stat says no file exists. Dying now."<<std::endl);
    } else {
        MMBLOG_FILE_FUNC_LINE(INFO, " stat says file exists. All is good. "<<std::endl);
    }
}



std::string   trim(const std::string& str, 
                 const std::string& whitespace// = " \t"
                 )   
{ 
    const auto strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1; 

    return str.substr(strBegin, strRange);
    //return std::string ("hello");
}

bool vectorCompare(String myString, vector<String> & comparisonStringVector) {
    if (comparisonStringVector.size() == 0) { //cout<<", returning TRUE "<<endl ;
        return true;} // If we are comparing to an empty vector, return true.  This is in case no partner chains have been specified, in which case any chain will pass.
    for (size_t i = 0; i < comparisonStringVector.size(); i++) {
        //cout<<__FILE__<<":"<<__LINE__<<" comparing "<<myString<< " to comparisonStringVector["<<i<<"] : "<<comparisonStringVector[i];
        if (comparisonStringVector[i].compare(myString ) == 0) { //cout<<", returning TRUE "<<endl ;  
            return true; } else {//cout<<", returning FALSE";
            }
        //cout<<endl;
    }
    return false; // If no String in comparisonStringVector is the same as myString
};

BondMobility::Mobility stringToBondMobility(String bondMobilityString) {
       String myBondMobilityString =   bondMobilityString;
       myBondMobilityString.toUpper();
       BondMobility::Mobility myBondMobility;
       // Remember  Free = 1, Torsion = 2, Rigid = 3
       if ((myBondMobilityString).compare("RIGID") == 0)              {
           MMBLOG_FILE_FUNC_LINE(INFO, " Detected Rigid : >"<<bondMobilityString<<"< or >"<<myBondMobilityString<<"< "<<std::endl);
           myBondMobility = SimTK::BondMobility::Rigid;
           //std::cout<<__FILE__<<":"<<__LINE__<<" returning myBondMobility = >"<<myBondMobility<<"< "<<std::endl;
       }
       else if ((myBondMobilityString).compare("TORSION") == 0)       {
           MMBLOG_FILE_FUNC_LINE(INFO, " Detected Torsion :"<<myBondMobilityString<<std::endl);
           myBondMobility = SimTK::BondMobility::Torsion;}
       else if ((myBondMobilityString).compare("DEFAULT") == 0)       {
           MMBLOG_FILE_FUNC_LINE(INFO, " Detected Default :"<<myBondMobilityString<<std::endl);
           myBondMobility = SimTK::BondMobility::Default;}
       else if ((myBondMobilityString).compare("FREE")  == 0)         {
           MMBLOG_FILE_FUNC_LINE(INFO, " Detected Free :"<<myBondMobilityString<<std::endl);
           myBondMobility = SimTK::BondMobility::Free ;}
       else {
           MMBLOG_FILE_FUNC_LINE(CRITICAL, " At this time only Default, Free,Torsion, and Rigid bondMobilities are supported. You are attempting to apply \"" << myBondMobilityString <<"\". "<<std::endl);
       }
       return myBondMobility;
}


void InterfaceContainer::addInterface(vector<String> myChains,vector<String> partnerChains,  double myDepth ,  String myMobilizerString ){
    Interface myInterface; 
    myInterface.Chains.clear(); myInterface.PartnerChains.clear();
    for (size_t i = 0; i < myChains.size(); i++) {myInterface.Chains.push_back( myChains[i]);}
    for (size_t i = 0; i < partnerChains.size(); i++) {myInterface.PartnerChains.push_back( partnerChains[i]);}
    myInterface.Depth = myDepth; 
    if (! (myMobilizerString.compare("NONE")==0) ||
	  (myMobilizerString.compare("Rigid")==0) ||
	  (myMobilizerString.compare("Default")==0) || 
	  (myMobilizerString.compare("Torsion")==0) ||
	  (myMobilizerString.compare("Free")==0) 
	) {
       MMBLOG_FILE_FUNC_LINE(CRITICAL, " Expected a mobilizer type (Default, Rigid, Torsion, Free), but got : >"<< myMobilizerString<<"< "<<std::endl);
    }
    myInterface.MobilizerString = myMobilizerString; 
    interfaceVector.push_back(myInterface); 
};


vector<TwoAtomClass> InterfaceContainer::retrieveCloseContactPairs(vector<MMBAtomInfo> & concatenatedAtomInfoVector ){
//vector<TwoAtomClass> InterfaceContainer::retrieveCloseContactPairs(  BiopolymerClassContainer & myBiopolymerClassContainer){
        OpenMM::NeighborList neighborList;
        openmmVecType boxSize = openmmVecType(10000,10000,10000);
        //vector<MMBAtomInfo> concatenatedAtomInfoVector = myBiopolymerClassContainer.getConcatenatedAtomInfoVector();
        vector<TwoAtomClass> contactingAtomInfoPairVector;
        contactingAtomInfoPairVector.clear();
        vector<openmmVecType> particleList(concatenatedAtomInfoVector.size());
        vector<set<int> > exclusions( particleList.size() );
        for (size_t i = 0; i < concatenatedAtomInfoVector.size() ; i++) {
            particleList[i] = concatenatedAtomInfoVector[i].position;
        }

        MMBLOG_FILE_FUNC_LINE(INFO, " neighborList size is : "<<neighborList.size()<<endl);
        for (int h = 0 ; h < numInterfaces(); h++ ){ // loop through interfaceContainer interfaces ..
            vector<String> referenceChains = getInterface(h).getChains();  
            vector<String> partnerChains = getInterface(h).getPartnerChains();  
            double         radius        = getInterface(h).getDepth();  
            MMBLOG_FILE_FUNC_LINE(INFO, "Now turning interface "<< h << " to individual constraints between pairs of atoms."<<endl);
            getInterface(h).print(); 
            cout<<__FILE__<<":"<<__LINE__<<endl;
            computeNeighborListVoxelHash(neighborList, particleList.size() , particleList, exclusions, &boxSize, false, radius  , 0.0);
            for ( size_t j = 0 ; j < neighborList.size(); j++) {
                if ((((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain , (referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(  partnerChains))) == 1))  || 
//Use an XOR here. This means if the 'partnerChains' evaluation is later set to return 1 when partnerChains is empty, this will still work.
                    ((( vectorCompare(concatenatedAtomInfoVector[neighborList[j].second].chain ,(referenceChains))) == 1) &&
                     (( vectorCompare(concatenatedAtomInfoVector[neighborList[j].first].chain  ,(  partnerChains))) == 1))     //Make sure that exactly one residue is in the 'referenceChains', and the other residue is in the 'partnerChains' .. thus only the desired interface is included
                                                                                                                         ) 
                     && (concatenatedAtomInfoVector[neighborList[j].first].chain.compare(concatenatedAtomInfoVector[neighborList[j].second].chain) != 0 )
                   ) // lastly,make sure the two atoms are not in the same chain.

                {   
                    TwoAtomClass myTwoAtomClass(
			(concatenatedAtomInfoVector[neighborList[j].first].chain),
			(concatenatedAtomInfoVector[neighborList[j].first].residueID),
			(concatenatedAtomInfoVector[neighborList[j].first].atomName),
			(concatenatedAtomInfoVector[neighborList[j].second].chain),
			(concatenatedAtomInfoVector[neighborList[j].second].residueID),
			(concatenatedAtomInfoVector[neighborList[j].second].atomName)//,
                        //(concatenatedAtomInfoVector[neighborList[j].first ].position - concatenatedAtomInfoVector[neighborList[j].second].position) // later, compute distance
                    );
                    contactingAtomInfoPairVector.push_back(myTwoAtomClass );
                    MMBLOG_FILE_FUNC_LINE(INFO, " Detected contact : ");
                    myTwoAtomClass.print();
                }
                else {
                    // Do nothing.
                }
            }
        }
        return contactingAtomInfoPairVector;  
};

ConstraintClass::ConstraintClass(){
        chain1 = ""; residueID1 = ResidueID(); atomName1 = ""; 
        chain2 = ""; residueID2 = ResidueID(); atomName2 = ""; 
        constraintType = WeldToGround ; 
        };  
ConstraintClass::ConstraintClass(String myChain, ResidueID inputResidueID,String myAtomName) {
        residueID1 = (inputResidueID);
        atomName1 = myAtomName;
        chain1 = myChain;
        residueID2 = ResidueID();
        atomName2 = "" ;    
        chain2 = ""; 
        constraintType = ( WeldToGround );
        //toGround = true ;
    }; 

ConstraintClass::ConstraintClass(String myChain, ResidueID inputResidueID,String myAtomName,String myChain2, ResidueID inputResidueID2,String myAtomName2, ConstraintType myConstraintType) {
        residueID1 = (inputResidueID);
        //residueID1.setInsertionCode ( residueID1.getInsertionCode());
        atomName1 = myAtomName;
        chain1 = myChain;
        residueID2 = (inputResidueID2);
        atomName2 = myAtomName2;
        chain2 = myChain2;
        setConstraintType(myConstraintType);
        //toGround = false;
    };
/*
Array_<MobilizedBodyIndex> ConstraintClass::fetchMobilizedBodyIndexArray_(BiopolymerClassContainer myBiopolymerClassContainer,SimbodyMatterSubsystem & matter ) {
        Array_< MobilizedBodyIndex >    coordMobod(2);
        coordMobod[0] =  myBiopolymerClassContainer.updBiopolymerClass(chain1).getAtomMobilizedBodyIndex(matter,residueID1    ,atomName1    );
        coordMobod[1] =  myBiopolymerClassContainer.updBiopolymerClass(chain2).getAtomMobilizedBodyIndex(matter,residueID2    ,atomName2    );
        return coordMobod;
    };
Array_<MobilizerQIndex> ConstraintClass::fetchMobilizerQIndexArray_(BiopolymerClassContainer myBiopolymerClassContainer,SimbodyMatterSubsystem & matter, State & state) {
        Array_< MobilizerQIndex >       coordQIndex(2);
        coordQIndex[0] =  MobilizerQIndex(0); //mobilizedBody1.getMobilizerQIndex(state);
        coordQIndex[1] =  MobilizerQIndex(0); // // mobilizedBody2.getFirstQIndex(state);
        return coordQIndex;
    }
*/

void ConstraintClass::setConstraintType (ConstraintType myConstraintType) 
{
    constraintType = myConstraintType;
}

ConstraintType ConstraintClass::getConstraintType () const
{
    return constraintType ;
}

String ConstraintClass::constraintTypeString () const  
{  // const promises not to change the object, i.e. ConstraintClass
        if (constraintType == WeldToAtom) { return  "WeldToAtom" ;}
        else if (constraintType == WeldToGround) {  return  "WeldToGround" ;}
        else if (constraintType == CoupledCoordinate) {  return  "CoupledCoordinate" ;}
        else if (constraintType == Undefined) { return "Undefined" ;}
        else return " ERROR! ";
}


void  ConstraintClass::print() const {
        std::cout<<__FILE__<<":"<<__LINE__   // to here is fine
          <<" : Chain ID : "      <<getChain1()
          <<" Residue    ID: "    <<getResidueID1().outString()
          <<" atom name : "       <<getAtomName1()
          <<" : Chain ID2 : "     <<getChain2()
          <<" Residue    ID2: "   <<getResidueID2().outString()
          <<" atom name2 : "      <<getAtomName2()
          //<<" to Ground: "        <<getToGround()
          <<" constraintType : " << constraintTypeString()
          <<endl;
    };  




int IntLen(const char* cstr)
{
  int    k, n = 0;
  if (cstr)
  {
    n = strspn(cstr,spaces);
    cstr += n;
    if (*cstr == '-' || *cstr == '+')
        ++cstr, ++n;
    k = strspn(cstr,digits);
    n = k?n+k:0;
  }
  return n;
}
/// <int>::[spaces][+|-]<digits>[garbage]
bool isNumber(string      line)
{
  //std::string test("1234.56");
  std::istringstream inpStream((line));
  double inpValue = 0.0;
  if (inpStream >> inpValue)
  {
    std::cout<<__FILE__<<":"<<__LINE__<<" Decided that "<<line<<" IS a number."<<std::endl; return 1;
    // ... Success!!  test is a number.
  }
  else
  {
    std::cout<<__FILE__<<":"<<__LINE__<<" Decided that "<<line<<" is not a number."<<std::endl; return 0;
    // ... Failure!!  test is not a number.
  }
}


bool isFixed (const String putativeFixedFloat) { // This checks that the string represents a floating point number in fixed format .. no scientific notation or other stray characters.
	int dotCount = 0;
	for (int i = 0 ; i < putativeFixedFloat.length() ; i ++) {
                //std::cout<<__FILE__<<":"<<__LINE__<<" putativeFixedFloat[i] = >"<<putativeFixedFloat[i]<<"< "<<std::endl;
		if (!((String(putativeFixedFloat[i]).compare("0") == 0) || 
		   (String(putativeFixedFloat[i]).compare("1") == 0) || 
		   (String(putativeFixedFloat[i]).compare("2") == 0) || 
		   (String(putativeFixedFloat[i]).compare("3") == 0) || 
		   (String(putativeFixedFloat[i]).compare("4") == 0) || 
		   (String(putativeFixedFloat[i]).compare("5") == 0) || 
		   (String(putativeFixedFloat[i]).compare("6") == 0) || 
		   (String(putativeFixedFloat[i]).compare("7") == 0) || 
		   (String(putativeFixedFloat[i]).compare("8") == 0) || 
		   (String(putativeFixedFloat[i]).compare("9") == 0) || 
		   (String(putativeFixedFloat[i]).compare(".") == 0) || 
		   (String(putativeFixedFloat[i]).compare("+") == 0) || 
		   (String(putativeFixedFloat[i]).compare("-") == 0)))  {
	 	        MMBLOG_FILE_FUNC_LINE(CRITICAL, "Found a character : >"<<putativeFixedFloat[i]<< "< which is not of [0-9].+- .. this is not a fixed/float!"<<endl);
			return false; // actually we shouldn't get to this line
		}
		if (((String(putativeFixedFloat[i]).compare("+") == 0) || (String(putativeFixedFloat[i]).compare("-") == 0)) && (i > 0)) {
	 	        MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have tried to use '+' or '-' somewhere other than at the beginning of the number string.  This is not allowed!"<<endl);
			return false; // actually we shouldn't get to this line
		}
		if  (String(putativeFixedFloat[i]).compare(".") == 0) {
			dotCount++;
			if (dotCount > 1) {
			    MMBLOG_FILE_FUNC_LINE(CRITICAL, "You have tried to use more than one '.' .. This is not allowed!"<<endl);
			    return false; // actually we shouldn't get to this line
			}
		} 
	} // of for
	return true;
}
/*
// a recursive algorithm for reading a double from a String.  This String may contain ints, user variables (begin with @), +, and -.  No whitespaces or additional characters should be in the String.
double   myAtoF(  map<const String,double> myUserVariables,  const char* value){
    MMBLOG_FILE_FUNC_LINE(INFO, "inside myAtoF. converting string : >"<<value<<"<"<<endl);

#ifdef Lepton_USAGE
    map<string,double> leptonFormatUserVariables; // Wish this were not necessary. but userVariables uses the signature const SimTK::String,double . Lepton uses string, double.
    leptonFormatUserVariables.clear();
    for (auto  myUserVariablesIterator = myUserVariables.begin() ; myUserVariablesIterator !=myUserVariables.end(); myUserVariablesIterator++) {
        leptonFormatUserVariables[myUserVariablesIterator->first] = myUserVariablesIterator->second;
    }
    double leptonResult = Lepton::Parser::parse(std::string(value)).evaluate(leptonFormatUserVariables);
    MMBLOG_FILE_FUNC_LINE(INFO, " Lepton evaluation = >"<< leptonResult<<"< "<<std::endl);
    return leptonResult; // if Lepton_USAGE is defined, then we return here and the rest of the procedure is not used. 
#endif
    // If Lepton_USAGE is NOT defined, then we parse the formula the old dumb way, as follows.

    size_t plusPosition  = String(value).find_last_of('+'); // returns the position of the last '+' in value
    size_t minusPosition = String(value).find_last_of('-'); // ditto for '-'
    if ((plusPosition > minusPosition) && (plusPosition  != String::npos) )  minusPosition = String::npos; // If the plus sign is closer to the end of the string, pretend we didn't find any '-'
    if ((plusPosition < minusPosition) && (minusPosition != String::npos) )  plusPosition  = String::npos; // Conversely, if the '-' is closer to the end, pretend we didn't find any '+' .. actually we might not have found any '+' anyway.
    String baseDoubleString ;
    double          increment = -1111;
    double          decrement = -1111;
    MMBLOG_FILE_FUNC_LINE(INFO, endl);
    if (plusPosition != String::npos) { // We have a '+' to deal with
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        baseDoubleString = String(value).substr(0, (plusPosition + 0) );
        String incrementString = (String(value).substr(plusPosition+1,1000)); // the second parameter is ridiculously large, but will be truncated at the end of the input String.
        //MMBLOG_FILE_FUNC_LINE(" The increment String is : "<<incrementString<<endl;
        stringstream incrementStringStream(incrementString);
        increment = myAtoF(myUserVariables, incrementString.c_str() );
        decrement = 0;
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
    } else if (minusPosition != String::npos ){ // we have a '-' to deal with
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        if (minusPosition== 0) {
            // If this is just a leading '-' sign, then put a zero to the left of that minus.
            MMBLOG_FILE_FUNC_LINE(INFO, "Detected a leading \'-\' sign. Will insert a zero to the left of the \'-\'."<<endl);
            baseDoubleString = "0.0";
        }
        else {
            // Otherwise, parse whatever is to the left of the minus sign:
            baseDoubleString = String(value).substr(0, (minusPosition + 0) ); }
        MMBLOG_FILE_FUNC_LINE(INFO, "baseDoubleString =  >"<<baseDoubleString  <<"< "<<endl);

        String decrementString = (String(value).substr(minusPosition+1,1000)); // the second parameter is ridiculously large, but will be truncated at the end of the input String.
        stringstream decrementStringStream(decrementString);
        MMBLOG_FILE_FUNC_LINE(INFO, "About to extract numerical decrement from the string >"<<decrementString<<"< "<<endl);
        decrement = myAtoF(myUserVariables, decrementString.c_str() );
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        increment = 0;
        //MMBLOG_FILE_FUNC_LINE(endl;
    } else { // no + or - found. This means we can return a result without further recursion.. i.e. we are at a leaf of the recursion tree.
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
        if (!((increment == -1111 ) && (decrement == -1111 )  )) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, "Unexplained error!"<<endl);
        }
        baseDoubleString = String(value);
        increment = 0;
        decrement = 0;
        double baseDouble;
        {
            MMBLOG_FILE_FUNC_LINE(INFO, endl);
            if ((baseDoubleString.substr(0,1)).compare("@") ==0) {
                MMBLOG_FILE_FUNC_LINE(INFO, endl);
                if (myUserVariables.find(baseDoubleString.c_str()) == myUserVariables.end())
                {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, "Undefined user variable "<<value<<endl);
                }
                MMBLOG_FILE_FUNC_LINE(INFO, "Read user variable "<<baseDoubleString.c_str()<<"  which is set to : "<<myUserVariables[baseDoubleString.c_str()]<<endl);
                baseDouble = double(myUserVariables[baseDoubleString.c_str()]);
            }
            else {
                // This has no '+', '-', or leading '@'.  However it's still possible that the user gave scientific notation, so e.g. 1e-5 leads here to baseDoubleString = '1e'.  So we need to make sure there is nothing but [0-9],'+','-' in this string:
        if (isFixed(baseDoubleString)) {
                    baseDouble = (atof(baseDoubleString.c_str()));
                        MMBLOG_FILE_FUNC_LINE(INFO, "We appear to be in a leaf of the recursion tree for myAtoF. Parsed >"<<baseDoubleString<< "< as : "<<baseDouble<<endl);
        } else {
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, "There was an error processing a putative floating point number."<<endl);
        }
            }
        }
        //MMBLOG_FILE_FUNC_LINE(endl;
        return baseDouble;
    }

    //MMBLOG_FILE_FUNC_LINE(endl;
    double baseDouble = myAtoF(myUserVariables,baseDoubleString.c_str() ) ;

    double finalDouble = baseDouble + increment - decrement;
    MMBLOG_FILE_FUNC_LINE(INFO, "Result of >"<< value  <<"< is : " << finalDouble <<endl);
    return finalDouble;
}*/
/*
bool aToBool(  const char* value ) {

    String upperValue(value);
    for(int i=0;i<(int)upperValue.length();i++)  {
        upperValue[i] = toupper(value[i]);
    }

    if (( upperValue ==  "TRUE" ) ||( upperValue ==  "1")) {
        MMBLOG_FILE_FUNC_LINE(DEBUG, "TRUE"<<endl);
        return true;
    }
    else if (( upperValue ==  "FALSE" ) ||( upperValue ==  "0")){
        MMBLOG_FILE_FUNC_LINE(DEBUG, "FALSE"<<endl);
        return false;
    }
    else {
        MMBLOG_FILE_FUNC_LINE(INFO, "Error -- you have specified"<<value<<endl);
        SimTK_ERRCHK_ALWAYS((upperValue == "TRUE" || upperValue == "FALSE" || upperValue == "1"  || upperValue == "0") ,"[ParameterReader.cpp]"," requires either True or False but was set to something else");
        return false;
    }

}
bool aToBool( const String& name, const char* value ) {
    return aToBool(value);
}
bool compareUpper( const String& param, const char* symbol ) {

    String upperParam(param);
    String upperSym(symbol);

    if( upperParam.length() != upperSym.length() ) return false;


    for(int i=0;i<(int)upperParam.length();i++)  {
        upperParam[i] = toupper(param[i]);
        upperSym[i] = toupper(symbol[i]);
    }

    if( upperParam ==  upperSym )
        return true;
    else
        return false;
}
*/


    // a recursive algorithm for reading an integer from a String.  This String may contain ints, user variables (begin with @), +, and -.  No whitespaces or additional characters should be in the String.
    int   myAtoI(  map<const String,double> myUserVariables,  const char* value){
        MMBLOG_FILE_FUNC_LINE(DEBUG, ""<<endl);
        size_t plusPosition  = String(value).find_last_of('+');
        size_t minusPosition = String(value).find_last_of('-');
        if ((plusPosition > minusPosition) && (plusPosition  != String::npos) )  minusPosition = String::npos;
        if ((plusPosition < minusPosition) && (minusPosition != String::npos) )  plusPosition  = String::npos;
        String baseIntegerString ;
        int          increment = -1111;
        //int          decrement = -1111;
        size_t lastPlusOrMinus = min(plusPosition,minusPosition);
        if ((lastPlusOrMinus != String::npos) && (lastPlusOrMinus != 0)) {
            baseIntegerString = String(value).substr(0, (lastPlusOrMinus + 0) ); // Put everything to the left of the last +/- into baseIntegerString
            String incrementString = (String(value).substr(lastPlusOrMinus+0,1000)); //  NOT adding 1 to lastPlusOrMinus means that the sign at lastPlusOrMinus goes with incrementString.
            MMBLOG_FILE_FUNC_LINE(DEBUG, " At this stage, we are adding >"<<baseIntegerString<<"< and >"<<incrementString<<"<"<<endl);
            increment = myAtoI(myUserVariables, incrementString.c_str() );
            MMBLOG_FILE_FUNC_LINE(DEBUG, incrementString<<" was interpreted as >"<<increment<<"<"<<endl);
        } else if (lastPlusOrMinus == 0) { // There is a leading + or -, and  this is the only +/- in the whole expression.
            baseIntegerString = String(value).substr(1, 1000); // Put everything to from position 1 onwards into baseIntegerString
            if (plusPosition == 0) {
                MMBLOG_FILE_FUNC_LINE(DEBUG, " Detected that the base string : >"<<String(value) <<"< has a leading \'+\'. "<<endl);
                // Trim the leading '+' and return the rest    
                return myAtoI(myUserVariables, baseIntegerString.c_str() );
                //cout<<__FILE__<<":"<<__LINE__<<" Interpreted >"<<baseIntegerString<<"< as "<<increment<<endl; 
	    } else if (minusPosition == 0) {
                MMBLOG_FILE_FUNC_LINE(DEBUG, " Detected that the base integer string : >"<<baseIntegerString<<"< has a leading \'-\'.  Inverting sign."<<endl);
                // Trim the leading '-' and return the negative
                return -myAtoI(myUserVariables, baseIntegerString.c_str() );
                //cout<<__FILE__<<":"<<__LINE__<<" Interpreted >"<<baseIntegerString<<"< as "<<increment<<endl; 
            }
        }
        else { // no + or - found.
            MMBLOG_FILE_FUNC_LINE(DEBUG, ""<<endl);
            if (!((increment == -1111 ))){// && (decrement == -1111 )  )) {
                MMBLOG_FILE_FUNC_LINE(CRITICAL, " Unexplained error!"<<endl);
            }
            MMBLOG_FILE_FUNC_LINE(DEBUG, ""<<endl);
            baseIntegerString = String(value);
            increment = 0;
            //decrement = 0;
            int baseInteger;
                if ((baseIntegerString.substr(0,1)).compare("@") ==0) {
                    cout<<__FILE__<<":"<<__LINE__<<""<<endl;
                    if (myUserVariables.find(baseIntegerString.c_str()) == myUserVariables.end())
                    {
                        MMBLOG_FILE_FUNC_LINE(CRITICAL, ": Undefined user variable "<<value<<endl);
                    }
                    MMBLOG_FILE_FUNC_LINE(DEBUG, ""<<endl);
                    double  intCast   = double(int(myUserVariables[baseIntegerString.c_str()]));
                    MMBLOG_FILE_FUNC_LINE(DEBUG, ""<<endl);
                    double  doubleCast = double(myUserVariables[baseIntegerString.c_str()]);
                    MMBLOG_FILE_FUNC_LINE(DEBUG, ""<<endl);
                    MMBLOG_FILE_FUNC_LINE(DEBUG, " Read user variable "<<baseIntegerString.c_str()<<"  which is set to : "<<myUserVariables[baseIntegerString.c_str()]<<endl);
                    SimTK_ERRCHK_ALWAYS(( (intCast) == doubleCast  ) ,"[ParameterReader.cpp]","Expected an int and got a non-integer");
                    MMBLOG_FILE_FUNC_LINE(DEBUG, ""<<endl);
                    baseInteger = int(myUserVariables[baseIntegerString.c_str()]);
                    MMBLOG_FILE_FUNC_LINE(DEBUG, ""<<endl);
                }
                else if (isNumber(baseIntegerString.c_str()))
                {
                    double  intCast   = double(int(atof(baseIntegerString.c_str())));
                    double  doubleCast = double(atof(baseIntegerString.c_str()));
                    SimTK_ERRCHK_ALWAYS(( (intCast) == doubleCast  ) ,"[ParameterReader.cpp]","Expected an int and got a non-integer");
                    baseInteger = (atoi(baseIntegerString.c_str()));
                } else {
                    MMBLOG_FILE_FUNC_LINE(CRITICAL, " : What you have entered: >"<<baseIntegerString<<"< is neither a variable (starting with @) nor an explicit number."<<endl);
                }
            //cout<<__FILE__<<":"<<__LINE__<<" : Result of "<<value<<" is : " <<  baseInteger <<endl;
            return baseInteger;
        }

        int baseInteger = myAtoI(myUserVariables,baseIntegerString.c_str() ) ; 

        int finalInteger = baseInteger + increment ;//- decrement;
        MMBLOG_FILE_FUNC_LINE(INFO, " : Result of "<< value  <<" is : " << finalInteger <<endl);
        return finalInteger;
    }   




/*void printBiopolymerSequenceInfo(const Biopolymer & myBiopolymer) {
    for (int i = 0; i < myBiopolymer.getNumResidues(); i++) {
        cout<<__FILE__<<":"<<__LINE__<<" Residue type, number, and insertion code: "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getOneLetterCode() <<", "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getPdbResidueNumber()<<", "<<myBiopolymer.getResidue(ResidueInfo::Index(i)).getPdbInsertionCode()<<endl;
    }    
};*/


         String intToString(int i) {
		/*
                // if boost doesn't work, go back to the stringstream method:
                String s ; 
                std::stringstream out;
                out << i;
                s = out.str();
                return s;
		*/
		//return  boost::lexical_cast<String>(i);      
                return SimTK::String(i);
            };  

Vec3 ValidateVec3(Vec3 myVec3){

    if (! myVec3.isFinite()){
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " This Vec3 vector is infinite or not a number : "<<myVec3<<endl);
    }
    return myVec3;
}



int ValidateInt (const int myInt) {
    /* WARNING: There used to be some validation code which was
     * apparently commented out at some point. Review git history
     * if you want to see the original (disabled) code.
     */
    return myInt;
}



int ValidateNonNegativeInt (const int myInt) {
	
        if (!(myInt >= 0)) {
			MMBLOG_FILE_FUNC_LINE(CRITICAL, " Expected a nonnegative integer. "<<endl);
		}
		else return myInt;
}




double ValidateNonNegativeDouble(const double myDouble) {
	
    if (std::isnan(myDouble)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Not a number! "<<endl);
    }
    if (std::isinf(myDouble)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Not a number! "<<endl);
    }
    if (!(myDouble>= 0)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Expected a nonnegative Double   . "<<endl);
    }
    return myDouble;
}

double ValidateDouble(const double myDouble) {
	
    if (std::isnan(myDouble) || std::isinf(myDouble)) {
        MMBLOG_FILE_FUNC_LINE(CRITICAL, " Not a number! "<<endl);
    }
    return myDouble;
}


Real DotProduct(const Vec3 vecA, const Vec3 vecB) {
	ValidateVec3(vecA);
	ValidateVec3(vecB);
	return (vecA[0]*vecB[0] + vecA[1]*vecB[1] + vecA[2]*vecB[2] );
}


vector<String> readAndParseLine   (ifstream & inFile) {
        stringstream u;
	String inString;
 	getline (inFile,inString);
	//u.str(inString);
	//cout<<__FILE__<<":"<<__FUNCTION__<<":"<<__LINE__<<" read: "<<inString<<endl;
	vector<String> mystring;
	
        istringstream iss(inString);

    	do
    	{
        	String sub;
        	iss >> sub;
        	//cout << "Substring: " << sub << endl;
		mystring.push_back(sub);
    	} while (iss);

	return mystring;
}

String removeAllWhite (String &str)
{
    String temp;
    for (int i = 0; i < str.length(); i++)
        if (std::string(str)[i] != ' ') temp += std::string(str)[i];
    str = temp;
    return str;
}

vector<String> readAndParseOnColWidth   (ifstream & inFile, int columnWidth) {
        stringstream u;
	String inString;
 	getline (inFile,inString);
	//u.str(inString);
	//cout<<__FILE__<<":"<<__LINE__<<" read: "<<inString<<endl;
	vector<String> mystring;
	
        istringstream iss(inString);
    	//for (int j = 0; j < numColumns; j++)
        	String sub;
	int j = 0;
	do
    	{
		sub = inString.substr(j * columnWidth,columnWidth);
		removeAllWhite(sub);
        	//cout << "Substring: >" << sub <<"<"<< endl;
		if (sub.length()>0)  mystring.push_back(sub);			
		j++;
    	} while (sub.length()>0);

	return mystring;
}


    /*ParameterStringClass::ParameterStringClass( const String & paramsLine ){
        char * params = strdup( paramsLine.c_str() );
        char * token = strtok( params, " " );

        clear();
        while( token ){
            add( token );
            token = strtok( NULL, " " );
            std::cout<<__FILE__<<":"<<__LINE__<<" Added token : >"<<token<<"< "<<std::endl;
        }
        free( params );
    }*/

    void ParameterStringClass::validateNumFields(int correctNumFields) const{ // make sure we have the right number of parameters
        std::cout<<__FILE__<<":"<<__LINE__<<" This line contains "<<size()<< " elements. comparing to "<<correctNumFields<<"."<<endl;
        if ( size() < correctNumFields ){
            MMBLOG_FILE_FUNC_LINE(CRITICAL, ": You have not specified enough parameters for this command."<<endl);
        } else if ( size() > correctNumFields ) {
            MMBLOG_FILE_FUNC_LINE(CRITICAL, ": You have specified too many parameters for this command."<<endl);
        };  
    };  
    void ParameterStringClass::print() const {
        MMBLOG_FILE_FUNC_LINE(DEBUG, endl);
        for (int i = 0 ; i < size(); i++){
            MMBLOG_FILE_FUNC_LINE(DEBUG, stringVector[i]<<" ");
            //std::cout<<__FILE__<<":"<<__LINE__<<" "<<i<<" >"<<stringVector[i]<<"< "<<std::endl;
        };
        MMBLOG_FILE_FUNC_LINE(INFO, endl);
    };

    String ParameterStringClass::getString() const {
        std::stringstream ss;
        for (int i = 0 ; i < size(); i++){
            ss <<" "<<stringVector[i];
        };

        return ss.str();
    }

