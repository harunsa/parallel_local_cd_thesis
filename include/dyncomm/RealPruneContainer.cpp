#include "parlay/parallel.h"
#include "../adj_listweight/adj_listweight.h"
#include <fstream>
#include <cinttypes>
#include "PruneRunner.cpp"
#include <iomanip>

using namespace std;

class RealPruneContainer {


public:
    static string path;
    static bool FAST_ESTIMATE;
    static int64_t HASH_ATTEMPTS;
    static int64_t INTERVALS_TO_CHECK;
    static bool DISPLAY_TIMINGS;
    static bool PRINT_EIGS;
    Eigen::MatrixXd bounds;
    double estimate;
    int64_t T;
    vector<bool> sizePruned;
    vector<bool> sizeUpdated;
    vector<bool> timePruned;
    vector<bool> timeUpdated;
    //onbasedtime =false ; pr_estimate 0.016; detailed_timing = true ; pr_usegroups  = true ; pr_compute_all_eigs = false
    RealPruneContainer(const AdjListWeight& g, int64_t T, const string& path, const string& filename, bool oneBasedTime,
                       double est, bool useGroups, bool useFull);
    bool improveEstimate(double e);
    void updatePruning();
    void printPruningStats();
    bool isPrunedInterval(int64_t ts, int64_t range);
    bool isPrunedDuration(int64_t range);
    static double getFastEstimate(const AdjListWeight& g,int T, MatrixXd &bounds);
    RealPruneContainer();
    bool isPruned(int64_t ts, int64_t range);
};

bool RealPruneContainer::FAST_ESTIMATE = false;
int64_t RealPruneContainer::HASH_ATTEMPTS = 100;
int64_t RealPruneContainer::INTERVALS_TO_CHECK = 10;
bool RealPruneContainer::DISPLAY_TIMINGS = true;
bool RealPruneContainer::PRINT_EIGS = true;
//here used static {} mb change later

RealPruneContainer::RealPruneContainer() =default;

//g is TGraph in sourcecode
RealPruneContainer::RealPruneContainer(const AdjListWeight &g, int64_t T, const string &path1,
                                       const string &filename,
                                       bool oneBasedTime, double est, bool useGroups,
                                       bool useFull) {



    cout << "T =" << T;
    this->T=T;
    this->estimate=est;
    this->bounds= Eigen::MatrixXd (T,T);

    if(false) {
        PruneRunner pr = {};


            cout << "\t@@@@@@ Load prune runner:"<< "s\n";


        Prune p = pr.loadFile(path1 + filename);


            cout << "\t@@@@@@ Done loading:"<< "s\n";

        pr.precompute(p, this->bounds, false);

        if (FAST_ESTIMATE) {
            //get.FastEstimate, false for now, implement later
        }
        if (useGroups) {
            this->bounds = pr.getNewMixedBounds(p, 0.016, this->bounds, useGroups);
        }
        if (DISPLAY_TIMINGS) {
            //...
        }
    }else{
        MatrixXd a {
                {0.11096587616411985, 0.11786616340223947, 0.094832956945698, 0.11850755377521324, 0.11861875265996995, 0.10401671874968721, 0.09802292348655661, 0.12039675208787413, 0.12132090184372091, 0.030382927344080685, 0.030382927344080685, 0.030382927344080685, 0.030382927344080685, 0.044027888014348406, 0.03881251594903216, 0.022123382978627436, 0.022123382978627436, 0.022123382978627436, 0.022123382978627436, 0.02633885405877963, 0.027551816127998643, 0.017953073348767296, 0.017953073348767296, 0.017953073348767296, 0.017953073348767296, 0.017953073348767296, 0.03518010275587115, 0.036346213549734385, 0.07734746482235828, 0.0385938647762527},
                {0.12249915542226372, 0.08596459475158208, 0.09280294765556259, 0.08673705608201017, 0.09005174456152917, 0.08707112341705504, 0.09757388700625064, 0.05454810230806666, 0.05594019660525684, 0.01789682007199237, 0.01789682007199237, 0.01789682007199237, 0.01789682007199237, 0.020432147923396315, 0.021235380533603607, 0.0185597576379762, 0.016752143861202256, 0.015961967020240598, 0.015303263174635299, 0.01677159496312585, 0.01768926375568962, 0.01750887093012862, 0.017742568143487136, 0.018649288976967965, 0.019180282106872688, 0.022874201314918054, 0.024694104152138758, 0.025124735512817437, 0.030345845545068844, 0.0},
                {0.12014220104823284, 0.11765796751566963, 0.09362216955694448, 0.09422890230305259, 0.08915104950604584, 0.1012617480745644, 0.09858411857358912, 0.03842703890455483, 0.029653041796165545, 0.03981338425726562, 0.030912836937640086, 0.025529440610946846, 0.022388195392987136, 0.01908850514553814, 0.017919462476612143, 0.015984362889464983, 0.015261279281668557, 0.014681455970483204, 0.015311228652684177, 0.016180876090584037, 0.01689277283696666, 0.017030777002732256, 0.017947522049683425, 0.01847955700026699, 0.022189187126408008, 0.02395777712458321, 0.024926796305945347, 0.03003043444212724, 0.0, 0.0},
                {0.11266873610148287, 0.07372933546419319, 0.08299216899739487, 0.08126793828836726, 0.0986992679371808, 0.09518864020178257, 0.09919162383608701, 0.026702980990053273, 0.022019647206348435, 0.02823617718983041, 0.023276717837272642, 0.020475425863260746, 0.018044099865206532, 0.017060710689308173, 0.01514825464275717, 0.01448359574054799, 0.013888840486768117, 0.014508481415391272, 0.01532616782301477, 0.016027683299225957, 0.01694035395613916, 0.017160568700520826, 0.017160568700520826, 0.017160568700520826, 0.017160568700520826, 0.017160568700520826, 0.029276364252513915, 0.0, 0.0, 0.0},
                {0.11648634733162402, 0.11857724785760013, 0.09418587189828831, 0.12162475410869428, 0.10562995147398704, 0.10875660246475045, 0.051467350664026836, 0.019741239259271515, 0.017067682250125622, 0.021246909901165886, 0.018765098282261423, 0.017270104213987603, 0.01637013316271284, 0.014387015743102267, 0.0137740072117511, 0.013294003925034693, 0.013888955534111716, 0.014651318794073881, 0.015387472793250496, 0.01625521261865577, 0.017188552114677732, 0.017775124905384698, 0.017775124905384698, 0.017775124905384698, 0.017775124905384698, 0.017775124905384698, 0.0, 0.0, 0.0, 0.0},
                {0.11837007081410064, 0.085821933087708, 0.05452128673619069, 0.09044955478684244, 0.12405727778688784, 0.04856505619640229, 0.030602035151268834, 0.01498182517094002, 0.013516668384635465, 0.016998812404990207, 0.016350222749026568, 0.015375803061780484, 0.013286527331973956, 0.012674744075472416, 0.012339174526805628, 0.012869500147017302, 0.013644015853519884, 0.014336819151888865, 0.014434369206999347, 0.016116878204246683, 0.01665652443767138, 0.017516980087391245, 0.017516980087391245, 0.017516980087391245, 0.017516980087391245, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.12068765882617169, 0.12317221614951915, 0.030217252091110168, 0.09919874479098037, 0.041501262147346174, 0.027441245120604002, 0.022775239949700116, 0.015072705380751802, 0.011646246729206569, 0.016028125706704498, 0.014881391973914236, 0.012713867286524964, 0.012149884482658798, 0.011855749222843099, 0.012342723798715718, 0.013005066195038752, 0.013688228976395931, 0.013735769652984606, 0.014661551959460343, 0.01611293077816316, 0.019924489932280954, 0.021873640666963526, 0.022603811096078775, 0.026895691409788234, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.12277015599982001, 0.08677468133343716, 0.01845591351108473, 0.033960036539625066, 0.025375412471814636, 0.021043024109540894, 0.018639978507984566, 0.009831709456440088, 0.015375416228550261, 0.014472766281080622, 0.011926072604545413, 0.011407551306691206, 0.011105467933157656, 0.011593092206118543, 0.012298188901053216, 0.012951961620629992, 0.012937989189451698, 0.013859204118447332, 0.014406916742593756, 0.019030202905067548, 0.021935660164500277, 0.02181864206384395, 0.025944261565709573, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.12652581733485546, 0.12535913814463528, 0.010486574124712633, 0.027224040571856678, 0.01997607702889685, 0.016564800583355997, 0.01172038128485334, 0.01610002839203785, 0.01448693354188049, 0.011554172144517606, 0.011075960427680937, 0.010745180928376563, 0.011186595813113607, 0.011787495179510478, 0.012403284160053025, 0.01229852752790319, 0.013207694857929164, 0.013745773564458903, 0.017519223934199576, 0.021140535188124022, 0.021871116882616636, 0.02511900496840371, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.12180316263398652, 0.01580901604866025, 0.006380502188431071, 0.010060571760722753, 0.009303637810714575, 0.00989506757363893, 0.00918839456866957, 0.006713496317338317, 0.011896657309351818, 0.013167805929542913, 0.012410899359122816, 0.012497430109707131, 0.013031839823794755, 0.011961215111772847, 0.01275522156961931, 0.014309673869177449, 0.014907601612765945, 0.01940964030016304, 0.021876038444028977, 0.020153836221493287, 0.023998256998838847, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.008351676359968682, 0.00917063363924685, 0.004183645949250249, 0.007332176032865695, 0.009192092814600874, 0.008006610419640995, 0.007418297419475871, 0.005920540682902903, 0.0058014097376092345, 0.007991502807974778, 0.008540574340543691, 0.009499583816490728, 0.010327077987172632, 0.01229578845712982, 0.013339763186153247, 0.016032446390888584, 0.016859219407201417, 0.020912614815484964, 0.02186917181919683, 0.023079319166121125, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.00989057161862586, 0.007187076761425574, 0.00442541304914695, 0.0070594396931596606, 0.008048425961627201, 0.007306642127447525, 0.007539081158704537, 0.0055390792314287726, 0.005981626823515799, 0.008043187029268598, 0.009161210884631944, 0.01018439133326424, 0.014513636606127023, 0.014792940811604333, 0.01586240092305102, 0.016391593111918023, 0.016391593111918023, 0.016391593111918023, 0.016391593111918023, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.00981819192872134, 0.009593295002380047, 0.0041136521682429835, 0.009340862836647333, 0.00766949866448638, 0.007822813357577266, 0.007311921842040534, 0.00589414749902877, 0.006927231591559637, 0.00946953374980426, 0.010641563343532152, 0.015702533844141275, 0.0172038099215577, 0.018498916117226514, 0.01940258749856834, 0.024487573975596227, 0.025626778253397063, 0.02674355789405732, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.009221310367896002, 0.005931942457034982, 0.003992376816559871, 0.006519899791385069, 0.006869785226405443, 0.006342937318597077, 0.007187762205211097, 0.006820939204038391, 0.008430628117156676, 0.011028167605487426, 0.016710690561516674, 0.01843610731202004, 0.019914627133239257, 0.02090249121067682, 0.02679106166060349, 0.02744100381095275, 0.02929682003368983, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.008208838199177284, 0.009016838142453224, 0.003880017384419421, 0.007364599170293308, 0.006646887410819835, 0.00766764757073723, 0.00880878816008759, 0.00893782976072955, 0.010623738121788339, 0.018602210836720468, 0.02056551611734216, 0.022155745628703993, 0.02351205901656105, 0.030361581565892137, 0.02991720897603508, 0.033153133068710997, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.009731112926194864, 0.0065051617211890666, 0.004024310270581134, 0.006549946497858549, 0.00904244850642303, 0.009598370103027372, 0.011531323412068854, 0.011832183234361762, 0.02104501138718094, 0.023207847054121027, 0.02378004298075787, 0.026785936808233813, 0.02838260159290065, 0.02995133982990284, 0.03834706174619161, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.009159065431842719, 0.008976073444624217, 0.0038519652988271203, 0.008840172151519369, 0.00973919856095732, 0.012194479589977863, 0.014500407964366914, 0.024924670178076494, 0.02524769759232859, 0.019728420016701522, 0.019728420016701522, 0.019728420016701522, 0.019728420016701522, 0.0387712517184056, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.008680891120777839, 0.0061787109730375885, 0.0038957788425315865, 0.00940719755421267, 0.013307600651822956, 0.016594407270628997, 0.02054764086445698, 0.022340797228947704, 0.02427066001624705, 0.02234672220189437, 0.02234672220189437, 0.02234672220189437, 0.02234672220189437, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.008418835985497114, 0.008634794424558466, 0.00866738440716297, 0.0170183585147886, 0.02149159904285545, 0.02727953750172181, 0.032737873509442204, 0.03129272059473453, 0.050993903436718904, 0.05755663478121831, 0.05725465522844122, 0.057255878132901734, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.00875590264353969, 0.015644203235440395, 0.020740129024705806, 0.033467852627909404, 0.041722809585719965, 0.04935632664404389, 0.05120350750868075, 0.04689677624052641, 0.0508166142170749, 0.06028920963924121, 0.06756078614797387, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.10514488380791989, 0.11215546280700657, 0.05173395409243395, 0.11301197835329345, 0.11727480139819668, 0.10271942353283497, 0.09805719085178687, 0.07687263029573149, 0.07584459102842882, 0.10251556353864129, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.11681369519027829, 0.08418949085549073, 0.05290418661908659, 0.09089116395773925, 0.09100052297775871, 0.050704589012060416, 0.09732263633144418, 0.08245478044247143, 0.09489151506453287, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.11291410981221965, 0.11297783018015302, 0.05261521786108135, 0.09342427932407966, 0.09166246450662066, 0.07565279083987547, 0.09570297213262638, 0.09820454616894263, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.11101663032102607, 0.08948204406909387, 0.05211848193145372, 0.08882260585605362, 0.10302628302984729, 0.0812276428684473, 0.09927827646676582, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.1301065635457264, 0.12349857232887389, 0.0557876365544785, 0.12267030281552575, 0.1031555377609815, 0.10687541255739476, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.11372996030218371, 0.08568407563144234, 0.05233601148733988, 0.08280842103657311, 0.11946815108307877, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.12027944791909706, 0.12030557971106759, 0.0658374546836663, 0.09284250404030275, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.11748649703842617, 0.07434600936744781, 0.08566140967049071, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.11088749227639151, 0.11844522571840349, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
                {0.12365945932082995, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
        };
        this->bounds=a;
    }



    sizePruned.resize(T);
    sizeUpdated.resize(T);
    timePruned.resize(T);
    timeUpdated.resize(T);

    //cout<<bounds<<endl;
}

bool RealPruneContainer::improveEstimate(double e) {

    if (e >= estimate) {
        return false;
    } else {
        estimate = e;

        for(int i = 0; i < T; ++i) {
            sizeUpdated[i] = false;
            timeUpdated[i] = false;
        }
        return true;
    }
}

void RealPruneContainer::updatePruning() {

    int t;
    bool chop;
    int size;
    for(t = 0; t < T; ++t) {
        sizeUpdated[t] = true;
        if (!sizePruned[t]) {
            chop = true;

            for(size = 0; size + t < T; ++size) {
                if (bounds(size,t) <= estimate) {
                    chop = false;
                    break;
                }
            }

            sizePruned[t] = chop;
        }
    }

    for(t = 0; t < T; ++t) {
        timeUpdated[t] = true;
        if (!timePruned[t]) {
            chop = true;

            for(size = 0; t + size < T; ++size) {
                if (bounds(t,size) <= estimate) {
                    chop = false;
                    break;
                }
            }
            timePruned[t] = chop;
        }
    }
}

bool RealPruneContainer::isPrunedInterval(int64_t ts, int64_t range){
    return bounds(ts,range) > estimate;
}

void RealPruneContainer::printPruningStats() {
    int total = 0;
    int pruned = 0;

    for(int t = 0; t < T; ++t) {
        for(int d = 0; t + d < T; ++d) {
            ++total;
            if (isPrunedInterval(t, d)) {
                ++pruned;
            }
        }
    }

    //4 decimals
    string perc = to_string(100.0 * ((double)pruned + 0.0) / ((double)total + 0.0));
    cout<<pruned<< " intervals pruned out of a total of "<<total<< " (" << perc << "%)."<<endl;

}

double RealPruneContainer::getFastEstimate(const AdjListWeight &g, int T, MatrixXd &bounds) {

    std::map<double, vector<int64_t>> eigs = {};
    double min = 1.0;
    int64_t st = true;
    int64_t dur = -1;

    for(int i = 0; i < T; ++i) {
        for(int j = 0; i + j < T; ++j) {
            if (bounds(i,j) > 0.0) {
                vector<int64_t> range = {};
                range.emplace_back(i);
                range.emplace_back(i + j);
                eigs.insert_or_assign(bounds(i,j), range);
                if (bounds(i,j) < min) {
                    min = bounds(i,j);
                    dur = j;
                }
            }
        }
    }

    double lamb = Prune::unnormalizeConductance(min, dur + 1);
    double upper = Prune::normalizeConductance(sqrt(4.0 * lamb), dur + 1);
    int attempt_counter = 0;

    if(HASH_ATTEMPTS>0){
        double bestCond = 1.0;

        for(auto it = eigs.begin();it!=eigs.end();++attempt_counter){

            double e = it->first;
            if (attempt_counter >= INTERVALS_TO_CHECK) {
                break;
            }

            vector<int64_t> alist = it->second;
            st = alist[0];
            dur  = alist[1] - alist[0];
            int64_t timebits = 2 * T / (dur+1);

            //hastable unordered set and map

            Edge a;//source , (end,weight) =  temporalNode
            ++it;
        }
    }
    return 0;
}

bool RealPruneContainer::isPrunedDuration(int64_t range) {
    if (sizeUpdated[range]) {
        return sizePruned[range];
    } else {
        sizeUpdated[range] = true;
        if (sizePruned[range]) {
            return true;
        } else {
            for(int t = 0; t + range < T; ++t) {
                if (bounds(t,range) <= estimate) {
                    return false;
                }
            }
            sizePruned[range] = true;
            return true;
        }
    }
}

bool RealPruneContainer::isPruned(int64_t ts, int64_t range) {
    int start = ts - range;
    if (start < 0) {
        start = 0;
    }

    for(int s = start; s <= ts && s + range < this->T; ++s) {
        if (this->bounds(s,range) <= this->estimate) {
            return false;
        }
    }
    return true;
}




