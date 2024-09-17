export
    SchubertQuadrics
function SchubertQuadrics(α,β,γ)
    Schuberts_Triangle = [[1,3,9,17,21,21,17,9,3,1],
                           [2,6,18,34,42,24,18,6,2],
                            [4,12,36,68,68,36,12,4],
                              [8,24,72,104,72,24,8],
                              [16,48,112,112,48,16],
                                  [32,80,128,80,32],
                                    [56,104,104,56],
                                        [80,104,80],
                                            [92,92],
                                               [92]]
    if α+β+γ != 9
        return("α+β+γ must equal 9")
    end
    enumerative_count=Schuberts_Triangle[β+1][α+1]
    @var p[1:4,1:α] l[1:6,1:β] h[1:4,1:γ]
    @var D x[1:4,1:4]
    #We represent quadrics by 4x4 symmetric matrices
    X = Matrix{Expression}(Symmetric(x))
    #Below we construct the second and third exterior power of X
    cols₂ =[[1,2],[1,3],[1,4],[2,3],[2,4],[3,4]]
    Λ₂X = [det(X[I,J]) for I in cols₂, J in cols₂]

    cols₃ = [[2,3,4],[1,3,4],[1,2,4],[1,2,3]]
    S=[[1,0,0,0] [0,-1,0,0] [0,0,1,0] [0,0,0,-1]]
    Λ₃X= S*[det(X[I,J]) for I in cols₃, J in cols₃]*S
    
    Point_Conditions=[p'*X*p for p in eachcol(p)]
    Line_Conditions=[l'*Λ₂X*l for l in eachcol(l)]
    Plane_Conditions=[h'*Λ₃X*h for h in eachcol(h)]
    
    Affine_Chart=sum(randn(Float64,10).*unique(vec(X)))-1
    Det_Value=det(X)-D
    
    params=vcat(vec(p),vec(l),vec(h))
    Equations=System(vcat(Point_Conditions,
                          Line_Conditions,
                          Plane_Conditions,
                          Affine_Chart,
                          Det_Value), 
                          parameters=params)
    E=EnumerativeProblem(Equations)
    return(E)
end

#=

Witness = Dict(
#Records of maximal reality

#110/112 (3,4,2)
(3,4,2)=>[[[1.0, 1.1452570794739443, 0.7275030922052121, 0.6910442726420025],[1.0, -0.2867207604845223, 1.8175635552718183, -0.9096803901468512],[1.0, -0.6521297186625569, 1.3190792836370007, 0.08501921278282337]],[[-2.1836987240486923, 0.7960763266995825, 0.16734663308159803, -2.859928783520865, -3.725476715059885, 1.1389686399335341],[-0.16451800237988204, 0.1651695290965864, 0.24942708460651036, -0.28674007182217537, -0.44364775696137687, 0.010675737136324515],[-0.20065189137810915, 0.5955773474383705, -0.7919498894439039, -0.8006726877051571, 0.7672768163376525, 0.8827225812177872],[-1.0978222454380382, 0.3381564737725297, -0.9615375018826063, -1.2544455215631474, 0.875851881411695, 0.8289333114427216]],[[-0.1830998286280456, 1.2697942894875682, 0.7152887480705059, 0.7302786760990851],[2.6656235228512952, -2.1911408306311926, 0.6683780262053292, -0.8331704066194641]]],
#96/104 (2,6,1)
(2,6,1)=>[[[1.0, 1.9999934443544602, 8.000002313236996, 6.999990896973108], [1.0, 1.0000152808614877, 8.999988454508744, 2.000017779922293]], [[-1.0000033788314642, 1.0948708106927932e-5, 3.2372596625350973e-6, -2.112220063301187e-6, -2.762761064401215e-5, 2.956478408795949e-10], [-3.2052703552099597e-6, -1.0000064320920476, -8.972394232652614e-6, -4.070104140681234e-6, -3.1175392494517966e-11, 1.6669376709353208e-6], [-5.670718710145783e-6, -1.1342940827102711e-5, -0.9999774980284154, 7.035716161882812e-11, 8.16990038690807e-6, 3.93514814411796e-6], [-0.9999940329853958, -1.1064538190910731e-5, 1.000005058894347, 7.137078823268617e-6, 0.999984178435409, 1.8201586670734888e-5], [-8.825581485562376e-6, -1.0000061748844384, 0.9999965913547298, -4.472299947560365e-6, 1.329780566025296e-5, 1.0000024456210113], [1.000011542278885, -0.999997305901678, -8.339884785233936e-6, -0.9999932736686208, -2.302749481073709e-5, 1.4687434557530988e-5]], [[-0.0603980952271978, 0.27025907222908246, 1.3625227181076238, 0.07545161522636186]]],
# 74/80 (3,5,1)
(3,5,1)=>[[[1.0, 1.845567723635855, 1.4523395133796857, -1.4190077878613492], [1.0, 0.7380218454571581, 1.3620805670593972, -0.7950215969142923], [1.0, -2.676218549810475, 1.1143273205193247, 1.105284551098701]], [[-0.21063422723986816, -0.1517768541063575, -2.3307791509504745, -0.30307152100558743, -1.0744527021447963, 2.5794275626930143], [-0.7853914542655948, -0.18763392884896204, 0.5994620625854895, -0.19760458584527948, 0.8542129778450424, 0.05325100530764219], [-1.522093714438261, 0.8888560393214976, -1.7192291826835182, 1.7813643511820458, -3.300403147179249, -0.08474531282452342], [-0.0340431512778693, 0.22947255010183185, 1.3805661580944557, 0.10129848913207473, 0.6777275495762244, -0.46031587927253925], [-1.0972434994233016, 0.3263233667508659, 0.6334146845196821, -0.867771931495417, 1.3984778151490442, -0.9168570820586187]], [[-1.155841807206707, -0.675958684239566, -0.5405091360406622, 0.2466154136761639]]],
# 84/104 (1,7,1)
(1,7,1)=>[[[1.0, -1.2174744687442838, 0.570925069524299, 1.7736559724240268]], [[-0.22898856720224642, -1.7504120312065092, 0.5606418090167855, -1.1733049108056308, -0.23456431348197138, -4.665682646021656], [-0.23042679527031618, 1.0191608218884434, -1.0254655150178846, -0.38526082943297174, 0.4816309323859052, -0.4156968025973156], [0.8395131994605569, -1.1461832056724346, 0.6042026577704246, 0.5475500183284282, -0.7207818293555207, 0.5900048405946986], [-0.9444463266273267, 0.6400979948963604, 0.7381906822266725, -0.6422280386303961, 3.0268609243508022, -2.553426589206648], [0.9668181361251813, -0.562526448641624, 0.6018604513421566, 0.5471045665756923, 0.6345726230925439, -0.7097968685138489], [1.093468131197935, -0.6190857206669091, -0.44638829578406214, 1.5077317894902589, 0.7582442534876344, 0.18621085342012458], [-0.4469958845397738, -1.5258246340972454, -0.5753805328109113, -0.7450215816142964, -0.3645399355338408, -0.28535631663862593]], [[-0.3714901579633248, -0.10073144296981888, -0.30934568617319447, -0.6269340799826092]]],
# 84/92 (1,8,0)
(1,8,0)=>[[[1.0, 0.9869561592509287, -0.8172726995919156, -0.9944084917530136]], [[-0.9994975714045787, -11.62653732043166, -101.37706341246671, -33.82075141644625, -393.06186487448485, -1141.873690193632], [-1.0002813626703262, 27.69009492131504, -575.0485137892168, -191.6816866148112, 5307.976763359258, -36741.77364554942], [-0.9996010978847363, -19.272888465827037, -278.58935235478685, -92.91373550231292, -1790.3412348448326, -8623.709552464086], [-1.000089340904571, -9.425639769885683, -66.62853184156069, -22.203548363291418, -209.28621948896938, -493.22261963006645], [-0.9995826613522095, -20.27782964267623, -308.4061296138084, -102.83373284264744, -2084.977850672565, -10568.682844757608], [-1.0006116204755569, 5.831112469674574, -25.499453770750193, -8.492157054016978, 49.53529979671828, -72.25634466097273], [-1.0001334314346109, 8.785345849257267, -57.8877923319539, -19.286732908570947, 169.47187180340129, -372.3529325546825], [-1.0003170821967444, 20.404993271954744, -312.27088043745795, -104.0689941074048, 2123.786541064606, -10834.69817543287]], []],
#Witnesses for full reality
(9,0,0)=>[[[1,181,-143,-22],[1,-192,98,9],[1,72,-181,79],[1,10,199,-154],[1,-98,-61,-186],[1,89,89,124],[1,6,106,21],[1,0,0,0],[0,0,0,1]],[],[]],
(8,0,1)=>[[[1,181,-143,-22],[1,-192,98,9],[1,72,-181,79],[1,10,199,-154],[1,-98,-61,-186],[1,89,89,124],[1,6,106,21],[1,0,0,0]],[],[[0,0,0,1]]],
(7,0,2)=>[[[1,181,-143,-22],[1,-192,98,9],[1,72,-181,79],[1,10,199,-154],[1,-98,-61,-186],[1,89,89,124],[1,6,106,21]],[],[[1,0,0,0],[0,0,0,1]]],
(8,1,0)=>[[[1,181,-143,-22],[1,-192,98,9],[1,72,-181,79],[1,10,199,-154],[1,-98,-61,-186],[1,89,89,124],[1,6,106,21],[1,0,0,0]],[[0,0,0,0,0,1]],[]],
(7,1,1)=>[[[1,181,-143,-22],[1,-192,98,9],[1,72,-181,79],[1,10,199,-154],[1,-98,-61,-186],[1,89,89,124],[1,6,106,21]],[[1,0,0,0,0,0]],[[1,0,0,0]]],
(7,2,0)=>[[[1,181,-143,-22],[1,-192,98,9],[1,72,-181,79],[1,10,199,-154],[1,-98,-61,-186],[1,89,89,124],[1,6,106,21]],[[0,0,0,0,0,1],[1,0,0,0,0,0]],[]],
(3,6,0)=>[[[1, 2, 8, 7],[1,1,9,2],[2,5,3,1]],[[249999//250000, 0//1, -9//1000000, -1//100000, 9//1000000, 1//500000], [1//1000000, 999997//1000000, -1//125000, 3//500000, -1//250000, -1//125000], [-1//250000, 1//1000000, 1000001//1000000, 1//500000, -3//500000, 1//500000], [9//1000000, 1//200000, -1//125000, 99999//100000, 3//1000000, -9//1000000], [-1//500000, 1//125000, -1//125000, 1//200000, 1000007//1000000, -3//500000], [3//500000, -3//1000000, 3//1000000, 1//200000, 1//500000, 500001//500000]],[]],
(0,3,6)=>[[],[[2.123330796888376, 0.6912623262448587, 0.005125493043124041, 0.6039201649069075, 3.1747229472471243, 1.0320911768927516], [0.6213408105365589, 0.5487156924218648, -0.13020946670429723, 0.46187922080252786, -0.42339468296071064, -0.27711403578551475], [-2.3335287720259816, 0.026985066083956344, -0.42566527318049735, -2.301418543949435, -0.9905014915612335, 0.4312621783296918]],[[-13.720944271893881, -6.213583833839737, -1.3352135068959994, 0.9790493993277065], [-16.46028874216132, 1.447676432516861, 6.528780350107965, 8.74160382880954], [-5.9012055592692025, -0.566454417957713, 1.7940886586873908, 2.444282131863796], [-5.122837352025064, -2.112087321793184, 0.19490357332449512, 0.5033372325639817], [-18.496443384968305, -0.3865945889961324, 5.1015433695336565, 1.089121070761917], [-3.223989537949238, -1.7548474029853731, -1.5535614469350802, -0.6561254193033143]]],
(0,4,5)=>[[],[[2.301224502047096, 1.2457430074529494, 0.3245850601187, -4.67705534880616, -2.708537826351831, -0.8065444133667952], [0.010482403764983963, -1.0324850537507348, -1.4631769187207873, 2.647435271579958, 3.756053445789105, -0.42002395460198844], [1.3451489465457875, -0.7662593858950055, -0.4155849336459346, -1.960890307948332, -0.13790650539091853, -0.5272608035007298], [-1.3973907906625673, -0.15090803106670375, -1.238011501260353, -0.9568598840634316, 0.18063779228312393, 0.8672328766101705]],[[-3.0441262399696725, 0.31806732894210515, -2.294936893925828, -0.8230920441839931], [4.905788384288908, -3.381868419238704, 0.8153311407009844, -0.19041867727428413], [0.9863102794726878, -1.3593220978187546, 1.4190224132576366, 1.497706573894195], [3.628706911757814, -0.5774073000877742, 0.3080918657728405, -1.7850285940096475], [5.159440785572494, -3.613879389840293, -1.3265941737428875, -2.917891710094424]]],
(0,5,4)=>[[],[[1.3720522516765028, 1.1953590768511186, -0.382588432232574, 1.3617507828226776, -3.779425272517464, -2.912994168268017], [-0.9381396071998337, -1.4107006680239118, -2.6814152175821038, -0.8604366725628496, 1.1250418343135835, 4.15107221239752], [1.5175463915750098, 0.6843684382569575, 1.2345224617728892, -3.344510867462202, 1.2547878761148226, 3.286628360316295], [1.4417068642785522, -0.6520006497749711, -0.6291939732116446, -3.5735330131765917, -1.4016923296599306, -0.9256674559210847], [-1.9684422708053002, 0.4648636891887938, -0.20786751875725706, 0.1646506196396465, -4.128707662672587, 0.957640865492334]],[[-0.5527278979928266, -0.6946431335427222, 1.2185015140151338, -2.6741352327557695], [-1.1672043717823914, -0.5768520925175487, -0.24787105348335317, -0.8042100681939295], [0.07749660900763225, 0.7193007113166627, -0.23311916436045185, 0.23329044655839365], [1.535923917121052, -1.123919809328802, -0.6360072376838071, -0.3653287457993123]]],
(0,9,0)=>[[],[[0.999675574783112, -14.998966827845981, 169.0006942081918, 55.97829492359802, -839.7866687629062, 3136.579269973095], [0.996931301160291, -7.0005512791063, 37.00092532936442, 12.004317736025147, -83.9520158234356, 143.98136301730386], [0.9990064865973081, 1.0009462086639918, 0.9999469096879015, 0.0032608106433450984, 0.0005463760957006978, -0.0027164432671876854], [0.999844617231605, 9.000386736014253, 61.00175909794575, 19.987087241275624, 179.92331594644193, 400.1941286537424], [1.0010955422682566, 16.998507997261907, 217.00063318684568, 71.90434735329477, 1223.279055315237, 5184.949568695906], [1.0015436518832592, 24.998190380968836, 469.0010364860221, 155.71983879505856, 3896.670992684024, 24339.38599177174], [1.0009373543482984, 32.99885099725548, 816.9991692203264, 271.74466326060934, 8972.292910729859, 73990.83710803364], [0.9992211228792658, 40.99845249471411, 1261.0005322558231, 420.29695019044897, 17226.873492689363, 176417.88448072952], [0.999708250403028, 48.998724971248826, 1801.0025623044608, 600.135336846313, 29404.10843720552, 360023.57968402316]],[]],
(3,3,3)=>[[[1, 439//922, -347//271, 67//343], [1, -211//484, 153//346, 257//254], [1, -575//404, 131//320, -37//42]],[[92//159, 92//293, -120//307, -77//256, -76//391, -96//311],[-107//114, -18//383, 109//116, -37//217, -45//307, -47//264], [-365//302, -45//368, 172//209, 74//245, 25//62, 87//353]],[[193//182, 75//397, -244//631, 195//272],[91//307, -17//122, -553//837, 70//309],[919//295, 103//36, 1199//371, 57//176]]],
(4,0,5)=>[[[1.0, -0.03011381881412807, 0.952750740758348, 1.0578474018904487], [1.0, 0.4104813671327917, 0.4314131961646348, 0.3726989281849856], [1.0, 0.9929922063591193, 0.33954131060089304, 0.4137176388288107], [1.0, -0.5287072560821856, -1.7465498792034324, -1.505715412092522]],[],[[-3.048405790682657, -4.131613946235795, -1.4863788879856066, -0.26813977361559643], [-0.20063400833881942, 0.25405713468061925, 0.7542408623066996, -0.6400846034231468], [0.3419063187356392, -2.3250700264467636, -0.42118254449579884, 1.5772128713457678], [-0.35529437779029266, 2.3201687332716645, -1.285766441174061, -1.53240748667453], [-1.5722043359079605, 0.30922505308818043, -0.37094507320554304, -0.4797190781194417]]],
(4,3,2)=>[[[1.0, 0.19937696065645488, -0.8962780856037743, 0.7175895363389434], [1.0, -0.9655521006635179, 0.6239343482127574, -1.1452532362767096], [1.0, -0.43619271084010247, -0.7151190279222437, -0.4695342913170511], [1.0, -0.15045907653091015, 0.5661516847110798, 0.3182746794219176]],[[-0.11116580410380773, -0.47956688967360545, -0.5345973186412596, -1.0500681902712747, -0.9059784994199584, 1.1414152837236349], [-0.1996300435402079, -0.8586040218323554, -0.7049993839657298, -0.15070772037441252, -0.7141968539033003, -2.5395147550345545], [1.2863292091079006, 0.8500201180604327, 1.6970519146160936, -0.11494042860553622, 4.144914778686013, 2.890644788000794]],[[-0.879071478348913, 3.539513703453217, 3.348609684700031, 1.8276455504174975], [-1.360566072217785, -0.025610827862781745, -1.1623619924247985, -1.0827144579379575]]],
(1,2,6)=>[[[1.0, 1.5010475056242596, -0.3715596667233419, 1.0231818577116525]],[[-3.7415553520224005, 8.879063721474829, -9.474797293030193, -24.053989045113827, 23.628504069804475, 4.839612222364558], [-1.9218154567522556, 2.666752249694196, -0.31300438039840317, -7.0121252414162445, 2.262493965204156, -1.9974264134857982]],[[1.5447485816523248, 7.201706807273422, -5.384493798790159, 8.75991222450175], [0.5173244424049176, 2.385127492943532, -0.0796432338975109, 0.751901635922988], [3895.2898824795807, -1891.15881062156, -1238.6188928326926, 622.3015473841734], [-7.1067122914529826, -10.525813052861592, 6.161202919148966, -2.0729801423123595], [0.12284529369566266, -0.49601167862022205, 2.016626200637681, 1.8494510735489866], [2.8400312570259105, -15.481071261253122, 132.67626333055685, 103.0118479662903]]],
(1,3,5)=>[[[1.0, -2.0072577919774544, 0.3054362941399763, -7.177573596870725]],[[-2.187918516082563, 0.22615001503949006, 0.4280180788069203, -2.47953011375105, -0.3319528488777387, -0.4507537948128554], [11.846474533368529, 0.324426640236103, -0.7489024390658061, 7.977562982118544, -4.776290535827455, 0.3735170722155313], [-1.401172975168338, -0.40323434683740345, 3.8945716755050532, -1.923346706320325, 1.71063244981366, -4.853666154397916]],[[3.24119033190556, -3.079612529288061, 0.9857682075669817, 0.5348091061076293], [-3.3841460456524524, -3.337626334240124, 0.1406586453108655, 1.514488587766626], [-3.0656589473911673, -3.01409992335484, 10.1177286956896, -3.7927210986781543], [-0.4121717441113782, 0.95759379769555, 0.5753254415004788, -2.8111095397801273], [6.5296205052713985, 3.682151603429516, -15.86665067380469, -16.509876127856863]]],
(4,4,1)=>[[[1.0, -0.41394462448313263, 1.015207777973526, -1.978868735337276], [1.0, -0.6603716079007259, 0.7631074953835579, 0.3573438661223989], [1.0, -1.6061345351343361, -0.22841324720048237, -1.2595794991399987], [1.0, -2.6773240888774503, -1.3906428980737158, -3.1347391776012548]],[[-1.0386220241957318, -0.538684687976004, -0.6002155095436091, -0.19919787894598892, -0.16143521928704652, 0.03138675564562178], [-0.05409113453721309, -1.357779218772501, 2.096432299191858, -0.3454830337188574, 0.5700221390141684, 0.9184947645811578], [-17.96164141142844, 15.678213388016676, 15.109779897132084, -65.02161496195492, -29.51629539346599, -28.933854157563072], [0.6575606347654597, 0.6574118359279365, -1.8727877018479053, 0.6195589834397592, -1.61917148935315, 0.14575073109680586]],[[2.658214205713778, 90.19569757291436, 5.066648182669239, -35.8484500437644]]],
(2,1,6)=>[[[1.0, 2.0080633885772796, -5.469420560214948, 1.3900551383484916], [1.0, -0.6160752268658393, 0.10274387605702896, 0.9147649086052676]],[[-1.9916816355835891, 2.3444778952962078, 1.986757201518537, 1.368584973516061, -1.172373245019917, 2.7452425690549065]],[[-4.690052142344988, 4.115747184286548, 0.7935147999654473, -0.3388045011132228], [-24.25769345963963, 15.008417049087635, 21.245772076512427, -1.3853046879805908], [0.41671654984468337, -8.32998431413697, 2.601221113188579, -1.389533459094394], [-8.145852278851907, 0.2612425352791199, -12.187984518461251, 9.652188858452881], [0.37634030960232684, -1.0011555842978241, -0.10183401012983928, -0.5058285754749113], [11.850336376870917, 8.519155954663766, -16.193015928070352, 11.16279583643565]]],
(2,2,5)=>[[[1.0, -13.955525744176708, -60.48413892212048, 41.259084458175], [1.0, -0.34206683817023603, -1.2787265620619674, -1.609512537650531]],[[50.63781153075665, -9.121704065580623, -26.08839510313362, 37.31352261354189, -17.680623563531114, 22.40869623253835], [3.3262242549694996, 0.3350693284288718, 0.013980628021706965, 2.107017333360015, -1.691866665593164, -0.17928738632505695]],[[-8.72343347964006, 12.488176942997109, -15.990513528481975, 7.621079570504993], [-2.2563687923334053, 3.2411889641594733, 3.4272460726666076, -2.3901547823627562], [5.689948648822002, -3.760037246604791, -1.6737105549030231, 0.3468739650915871], [27.777080390915703, -0.29267314583918186, -18.089280423606933, -15.344624255156115], [37.27888974049697, 1.7588960391145287, -3.921777538535308, -23.89059699798142]]],
(3,0,6)=>[[[1.0, 24.59498006663152, -3.967962022680723, -0.850435073703259], [1.0, -12.435082328072665, 3.085359704893075, -1.8503776657123523], [1.0, -2.1871582679875976, -2.591302723550152, -3.588515882630248]],[],[[194.71186607356435, -12.145044772185406, -52.04202034998901, -53.85960185233909], [-0.662891010763335, -0.40497407492097975, -10.29327027439726, 5.121672683189074], [6.1713083408373315, 5.118020205741839, -0.6051185066670663, 1.1493206558616194], [-10.201795721740927, -11.48400603197178, -9.728724346949038, -1.5616243058333763], [156.80664249009064, -502.08014454864144, 111.32960541225323, -453.0672200693393], [-0.48208107801681543, -0.06031962610419911, -1.5057388119557045, 2.1344113563053737]]],
(3,1,5)=>[[[1.0, 0.11196402549538087, 0.8689036478993275, -1.2268063088062071], [1.0, -1.281099856638499, -0.4119814931004374, -0.3030724818582338], [1.0, 0.706141542951976, -0.5859420843064734, -0.6055773857168295]],[[-2.904621300713732, -2.8973876321973604, 3.275141120301843, 5.613460612417988, -2.284519885083963, 4.050695357595362]],[[0.5606856455243108, -6.209470128736699, -0.5011891526682949, -0.4140101952628523], [1.2236220748793238, 0.861869454947733, 1.6700329239361702, 1.8399581282577127], [3.2339874852929547, 0.8483605730062205, 1.0122117929170864, 5.595245076526597], [18.731185506111263, 27.349585135057893, -53.32980040024048, 53.416241094848374], [-0.7890550685583911, 2.9454465771589526, 2.9229547754595746, -1.1299420257387174]]],
(4,2,3)=>[[[1.0, -1.5025525842553848, 1.5575579758469125, -2.3909201450743907], [1.0, -0.5612682807313427, -0.17802117189951083, -0.23100174705960583], [1.0, -2.616068601926085, -2.0916936469993868, -5.108611357075233], [1.0, -0.18152152036486727, 1.2369486621535197, -4.021697393032103]],[[-1.8565794983214357, 4.82988933321377, -2.1791541816827085, 8.28406655015502, -7.493353105314875, 9.77055277495897], [-2.6833432019298424, 5.640296243926144, -3.494375893579112, -0.049576113954022054, -3.6711579144012245, 7.781209559007854]],[[30.584337029318505, 49.19811818690047, -44.63253991506339, 48.99083945293469], [-1.2383749099732464, -3.7412651383086257, 3.2645282300054626, 1.425226234863439], [4.218129922100462, -2.4513238466288256, 5.243610286807471, 24.609822216495978]]],
(4,1,4)=>[[[1.0, -4.559167667836678, 2.7987083385472826, -6.61309348392271], [1.0, 0.38824376025840157, 0.15039733146688847, -0.30913432156307125], [1.0, 12.791612607999227, 3.54948358970739, -5.821745376245071], [1.0, -6.093489180027147, -2.0870558286310668, -2.216694961244847]],[[0.49097975626410906, -1.3097951301898396, -3.5850912747793773, 1.1766261619406413, 3.1657840624892652, 0.14620488452660346]],[[-0.07536145412862749, 0.5591597903667405, -0.7130498830546519, 1.4489319726370589], [1.8532268504090363, -10.397934280836823, -17.72368466697672, -15.832150144582855], [-2.4186745262148746, 2.04792364791863, -0.6292285367771797, -4.523337573640416], [-4.738954914397866, 6.485415475363126, 5.993590162041933, -4.270029006727771]]],
(2,7,0)=>[[[1.0, 0.958706792420631, -0.46647129783076746, -0.3148809397288257], [1.0, 1.027372295820429, -0.32117476914165777, -0.2114059407498258]],[[0.4008645388607157, -1.2201890467178078, -0.48202081010915593, -0.8472154538477042, -0.16686068674198878, -0.5108306603828363], [0.609967737245868, 0.34724457611839793, 0.14045396189571213, -0.6194425578954955, 0.03394216500295152, 0.16195839237052406], [0.09688663789358037, -0.7231088525421742, -0.5841673463576765, 0.6620139820963818, 0.6280628265054468, -0.6959766598557385], [0.5769226239332654, -0.35119017557302606, -2.1452646360764773, 0.08591326557213597, -1.5925272576367921, 1.2888844825438115], [-0.015318050033647786, -0.32894789102316113, 1.5902404448560379, 0.0515702091701625, -0.18144948583327603, 1.457209411834611], [-1.4655451759813563, -0.6129371166251921, 0.7496750692260776, 0.9072681026070772, 2.342609908706074, 1.4438509812783806], [2.024703739060808, -0.27842745786755413, -1.7447810908547936, 1.7199734331941852, 5.122309499259092, 0.7777856485071597]],[]],
(2,5,2)=>[[[1.0, -0.7612529775987225, -0.9041707921308169, -0.5979597434804667], [1.0, -0.5324491029095367, -0.4767302423300307, -0.22331741618655]],[[-0.42059217344660077, 0.38522161456295695, -0.39068249286472456, 0.8271215152916483, -1.214450372652219, 0.34401647723921125], [0.7211281613404645, 1.1713097436380542, -0.6503493219545826, -0.24771202703588707, 1.2388791143642874, 1.7888801161555783], [0.9229140951819109, 0.7872267841881713, 0.7125212652798897, 0.8693212032292349, 0.5500333107979628, -0.20197859169865118], [-2.235915819532634, 0.7438078982674584, -0.001750215906302979, 1.345046337876959, 2.6491563074571056, -0.8823308505144387], [-0.6345408459316946, -0.20393234102409802, 2.5149848362898735, 0.08550613185730532, 2.7719959251405086, 1.2297809484719764]],[[0.36447147005511726, -1.1020983627771812, -0.48265305503504924, -0.7083613647089662], [2.292597407006857, -3.1683493679025623, -2.76808973219021, -2.0731208014290154]]]
)



=#