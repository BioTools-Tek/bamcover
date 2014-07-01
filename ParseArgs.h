#ifndef PARSEARGS_H
#define PARSEARGS_H

#include <QString>
#include <QStringList>
#include <iostream>

using namespace std;

class ParseArgs
{
public:
    ParseArgs(int argc, char *argv[]);
    //PARSE ARGS
    int depth;
    QString bamfile;
    QString optregion;
    QString bedfile;
    bool keepchr, debug;
};

void usage(char * name){
    char * version = "0.6";

    cerr << endl;
    cerr << "Prints percentage pileup in a given region for a given BAM file" << endl;
    cerr << "usage: " << name << " <in.bam> <in.bed> [OPTIONS]             v" << version<< endl;
    cerr << "       " << name << " <in.bam> <region> [OPTIONS]" << endl;
    cerr << endl;
    cerr << "OPTIONS:" << endl;
    cerr << "  --depth N  specifies that all positions with a read depth or coverage" << endl;
    cerr << "             less than N are deemed bad'. Default is 1." << endl;
    cerr << endl;
    cerr << "  --DEBUG  Prints manual pileup" << endl;
    cerr << endl;
    cerr << "  --chr    Use this if region not giving any data (preserves 'chr')" << endl;
    cerr << endl;
    cerr << "To manually check a pileup region, try:" << endl;
    cerr << " e.g. " << name << " myfile.bam chr2:123456-123457 --depth 27 --DEBUG" << endl;
    cerr << endl;
    cerr << "WARNING: make sure that a .bai index file is present, and that SAMTOOLS is installed" << endl;
    exit(0);
}

ParseArgs::ParseArgs(int argc, char *argv[]){
    //PARSE ARGS
    optregion="";
    bamfile="";
    depth = 1;
    keepchr=false;
    debug=false;


    if (argc<3) usage(argv[0]);

    if (argc >=3 ){
        bamfile = argv[1];
        QString second = QString(argv[2]);

        if (second.contains(":") && !second.contains("bed")) { //region
            optregion = argv[2];
        } else {
            bedfile = argv[2];
        }

        //Opts
        bool errs=false;
        QString errorsT = "Cannot parse: ";
        QStringList okay_vars;

        for (int i=3; i < argc; i++){
            QString tmp = argv[i];
            if(tmp=="--chr") keepchr=true;
            else if(tmp=="--DEBUG") debug=true;
            else if(tmp=="--depth"){
                if (argc <= i+1){
                    cerr << "Give a depth value. Default is 0" << endl;
                    exit(-1);
                }
                depth = QString(argv[i+1]).toInt();
                okay_vars.append(argv[i+1]);
            }
            else {
                bool error_yes=true;
                for (int i=0; i< okay_vars.length(); i++){
                    if (tmp==okay_vars.at(i)){
                        error_yes=false;
                        break;
                    }
                }
                if (error_yes) errorsT.append(tmp+" ");
            }
        }
        if (errs) {
            cerr << errorsT.toUtf8().data() << endl;
            exit(-1);
        }
    }
}

#endif // PARSEARGS_H
