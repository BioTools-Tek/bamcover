#include "ParseArgs.h"

#include <QProcess>
#include <QFile>
#include <QTextStream>


class LineHolder{
public:
    QString chrom;
    uint position;
    int depth;

    LineHolder(QString chrom, uint position, int depth){
        this->chrom = chrom;
        this->position = position;
        this->depth = depth;
    }
};


QProcess *qp = 0;

QList<LineHolder*> performPileup(QString bamfile, QString chrom, uint pos1, uint pos2, bool debug)
{
    qp = new QProcess;

    QString command = "samtools mpileup -A -r ";
    command.append(chrom).append(':').append(QString::number(pos1)+'-'+QString::number(pos2))\
            .append(' ').append(bamfile);

    cerr << command.toUtf8().data() << endl;

    qp->start(command);
    qp->waitForFinished();

    QStringList res = QString(qp->readAllStandardOutput()).split('\n');

    qp->terminate();
    qp->close();
    qp->kill();

    QList<LineHolder *> results;

    int i=0;
    uint last_pos = pos1;
    while (i < res.length()){
        if (res[i].length() > 3){

            QStringList bar = res[i].split('\t');
            //Assume constant chromosome -- fair assumption, bed files dont span chroms
            uint current_pos = bar[1].toInt();

            if (current_pos > pos2) {
                cerr << "breaking early!" << endl;
                break;
            }
            //Fill blank regions before
            while (last_pos < current_pos ){
                LineHolder *lh = new LineHolder(chrom, last_pos++, 0);
                results.append(lh);
            }
            //Add current
            results.append(new LineHolder(chrom, current_pos, bar[3].toInt()));
            last_pos++;
        }
        i ++;
    }
    //Fill blank regions after
    bool add_ends=false;
    while (last_pos < pos2){
        add_ends=true;
        LineHolder *lh = new LineHolder(chrom, last_pos++, 0);
        results.append(lh);
    }
    if (add_ends) results.append(new LineHolder(chrom, last_pos, 0));




    if (debug){
        cout << "\nchrom\tpos\tdepth" << endl;
        cout << "======================" << endl;
        for (int l=0; l< results.length() ; l++)
                cout << results[l]->chrom.toUtf8().data() << '\t' << results[l]->position << '\t' << results[l]->depth << endl;
        exit(0);
    }
    return results;
}


void printBed(QString bamfile, QString chrom, uint pos1, uint pos2, int depth, bool keepchr, bool debug, QString name="")
{
    //At this stage chrom should be a single number
    if (keepchr) chrom = "chr"+chrom;

    //DO PILEUP
    QList<LineHolder*> piles = performPileup(bamfile, chrom, pos1, pos2, debug);



    int length = piles.length();
    int above_depth=0;

    int i=0;

    bool recording = false;
    uint start_region_pos = 0, end_region_pos = 0;

    QString last_print="";

    while (i < length){
        LineHolder *pr = piles[i];
        bool printing=false;

        if (pr->depth >= depth){
//            cerr << pr->chrom.toUtf8().data() << " " << pr->position << "  " << pr->depth << endl;

            if (!recording) {
                recording = true;
                start_region_pos = pr->position;
            }
            else {
                //while recording update end position
                end_region_pos = pr->position;
            }
            above_depth ++;
        }
        // Found bad depth, stop recording, print current region
        else if (recording) printing = true;

        if (i >= length -1)
            if (start_region_pos!=0) printing = true;

        if (printing){
            QString complete = chrom+'\t'+QString::number(start_region_pos)+'\t'+QString::number(end_region_pos)+'\t'+name;
            if (last_print!=complete) cout << complete.toUtf8().data() << endl;
            last_print = complete;
            recording = false;
        }
        i++;;
    }

    cerr << chrom.toUtf8().data() << ':' << pos1 << '-' << pos2
         << '\t' << name.toUtf8().data() << '\t' << ((100*above_depth)/length) << '%' << endl;
}


void performPileupBed(QString bamfile, QString bedfile, int depth, bool keepchr, bool debug)
{

    QFile inputFile(bedfile);
    if (inputFile.open(QIODevice::ReadOnly))
    {
        QTextStream in(&inputFile);
        while ( !in.atEnd() ){
            QString line = in.readLine().trimmed();
            if (line[0]=='#'){   // print all headers out
                cout << line.toUtf8().data() << endl;
            }
            else{
                if (line.length()>0){
                    QStringList tokens = line.split('\t');

                    QString chrom = tokens[0];
                    if (chrom.startsWith("chr")) chrom = chrom.split("chr").last();
                    uint pos1 = tokens[1].toUInt();
                    uint pos2 = tokens[2].toUInt();
                    QString name = "   ";
                    if (tokens.length()>3) name = tokens[3];


                    if (pos1>=pos2){
                        cerr << "Skipping " << pos1 << "  " << pos2 << " line. Bad ordering." << endl;
                        continue; // Skip ambiguous or incorrectly ordered beds
                    }


//                    cerr << "CHROM: " << chrom.toUtf8().data() << endl;
                    printBed(bamfile, chrom, pos1, pos2, depth, keepchr, debug, name);
                }
            }
        }
    }
    inputFile.close();
}


int main(int argc, char *argv[])
{
    ParseArgs *pa = new ParseArgs(argc, argv);

    if (pa->optregion!="") {
        QStringList chr_regs = pa->optregion.split(':');
        if (chr_regs[0].startsWith("chr")) chr_regs[0] = chr_regs[0].split("chr").last();
        QStringList reg1_reg2 = chr_regs.at(1).split('-');

        printBed(pa->bamfile, chr_regs[0], reg1_reg2[0].toUInt(), reg1_reg2[1].toUInt(),pa->depth, pa->keepchr, pa->debug);
    }

    else performPileupBed(pa->bamfile, pa->bedfile, pa->depth, pa->keepchr, pa->debug);

    return 0;
}
