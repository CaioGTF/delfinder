import os
import subprocess
import warnings
from Bio import SeqIO
from Bio import Align
from Bio import AlignIO
from Bio.Align.Applications import MafftCommandline

lista = []
unique_ids = set()


felidae_count = 0
hyaenidae_count = 0
herpestidae_count = 0
odobenidae_count = 0
ursidae_count = 0
phocidae_count = 0
otariidae_count = 0
canidae_count = 0
mustelidae_count = 0
physeteridae_count = 0
balaenidae_count = 0
phocoenidae_count = 0
rhinocerotidae_count = 0
balaenopteridae_count = 0
monodontidae_count = 0
delphinidae_count = 0
pontoporiidae_count = 0
lipotidae_count = 0
camelidae_count = 0
ziphiidae_count = 0
hippopotamidae_count = 0
vespertilionidae_count = 0
suidae_count = 0
pteropodidae_count = 0
equidae_count = 0 
cervidae_count = 0
tapiridae_count = 0
bovidae_count = 0
rhinolopidae_count = 0
moschidae_count = 0
indriidae_count = 0
leporidae_count = 0
cheirogaleidae_count = 0
cercopithecidae_count = 0
dasypodidae_count = 0
hominidae_count = 0
sciuridae_count = 0
callitrichidae_count = 0


felidae = ['Herpailurus yaguarondi', 'Prionailurus bengalensis', 'Otocolobus manul', 'Felis catus', 'Acinonyx jubatus', 'Prionailurus viverrinus', 'Puma yagouaroundi', 'Lynx lynx', 'Panthera onca', 'Lynx rufus', 'Panthera tigris', 'Panthera pardus', 'Panthera uncia', 'Leopardus geoffroyi', 'Panthera leo', 'Neofelis nebulosa', 'Lynx canadensis']
hyaenidae = ['Hyaena hyaena']
herpestidae = ['Suricata suricatta']
odobenidae = ['Odobenus rosmarus']
ursidae = ['Ailuropoda melanoleuca', 'Ursus arctos', 'Ursus americanus', 'Ursus maritimus', 'Ursus americanus']
phocidae = ['Mirounga leonina', 'Mirounga angustirostris', 'Neomonachus schauinslandi', 'Halichoerus grypus', 'Phoca vitulina', 'Leptonychotes weddellii']
otariidae = ['Arctocephalus gazella', 'Callorhinus ursinus', 'Eumetopias jubatus', 'Zalophus californianus ']
canidae = ['Nyctereutes procyonoides', 'Vulpes vulpes', 'Canis lupus', 'Canis aureus', 'Canis latrans', 'Vulpes lagopus', 'Alopex lagopus', 'Canis familiaris']
mustelidae = ['Martes flavigula', 'Pekania pennanti', 'Mustela erminea', 'Meles meles', 'Lutra lutra', 'Enhydra lutris', 'Neogale vison', 'Lontra canadensis', 'Mustela lutreola', 'Mustela nigripes', 'Mustela putorius', 'Neovison vison', 'Martes pennanti', 'Martes zibellina', 'Mustela sibirica', 'Mustela altaica', 'Mustela nivalis', 'Martes melampus', 'Martes zibellina', 'Martes martes']
physeteridae = ['Kogia breviceps', 'Physeter catodon']
balaenidae = ['Eubalaena glacialis']
phocoenidae = ['Neophocaena phocaenoides', 'Neophocaena asiaeorientalis', 'Phocoena phocoena', 'Phocoena sinus', 'Phocoenoides dalli']
rhinocerotidae = ['Ceratotherium simum', 'Diceros bicornis', 'Diceros bicornis', 'Dicerorhinus sumatrensis', 'Rhinoceros unicornis', 'Ceratotherium simum']
balaenopteridae = ['Balaenoptera acutorostrata', 'Balaenoptera musculus', 'Megaptera novaeangliae', 'Balaenoptera physalus', 'Balaenoptera ricei']
monodontidae = ['Monodon monoceros', 'Delphinapterus leucas']
delphinidae = ['Sotalia fluviatilis', 'Sousa chinensis', 'Tursiops aduncus', 'Delphinus capensis', 'Tursiops truncatus', 'Lagenorhynchus albirostris', 'Orcinus orca', 'Orcaella heinsohni', 'Stenella longirostris', 'Lagenorhynchus obliquidens', 'Globicephala melas', 'Cephalorhynchus commersonii', 'Delphinus delphis', 'Lagenorhynchus obscurus', 'Lagenodelphis hosei', 'Steno bredanensis', 'Globicephala macrorhynchus', 'Pseudorca crassidens', 'Feresa attenuata', 'Stenella attenuata', 'Lagenorhynchus acutus', 'Peponocephala electra', 'Stenella coeruleoalba', 'Grampus griseus']
pontoporiidae = ['Pontoporia blainvillei']
lipotidae = ['Lipotes vexillifer']
camelidae = ['Vicugna vicugna', 'Vicugna pacos', 'Lama glama', 'Vicugna vicugna', 'Lama guanicoe', 'Camelus dromedarius', 'Camelus bactrianus', 'Mesoplodon bidens', 'Lama guanicoe', 'Lama pacos', 'Camelus ferus']
ziphiidae = ['Mesoplodon densirostris', 'Ziphius cavirostris', 'Mesoplodon bidens']
hippopotamidae = ['Hippopotamus amphibius']
vespertilionidae = ['Myotis davidii', 'Myotis daubentonii', 'Myotis lucifugus', 'Myotis myotis', 'Myotis brandtii', 'Eptesicus fuscus', 'Miniopterus natalensis']
suidae = ['Sus scrofa', 'Phacochoerus africanus', 'Potamochoerus porcus', 'Babyrousa babyrussa']
pteropodidae = ['Rousettus aegyptiacus', 'Pteropus vampyrus', 'Pteropus giganteus', 'Pteropus alecto']
equidae = ['Equus asinus', 'Equus caballus', 'Equus quagga', 'Equus przewalskii', 'Equus grevyi', 'Equus zebra', 'Equus hemionus', 'Equus kiang', 'Equus quagga', 'Equus burchellii']
cervidae = ['Odocoileus virginianus', 'Rangifer tarandus', 'Cervus canadensis', 'Odocoileus virginianus', 'Cervus elaphus', 'Capreolus capreolus','Muntiacus reevesi', 'Dama dama', 'Cervus canadensis', 'Cervus dama']
tapiridae = ['Tapirus bairdii', 'Tapirus terrestris', 'Tapirus pinchaque']
bovidae = ['Rupicapra pyrenaica', 'Ovibos moschatus', 'Rupicapra rupicapra', 'Budorcas taxicolor', 'Capra hircus', 'Alces alces', 'Oryx dammah', 'Ovis aries', 'Antidorcas marsupialis', 'Ovis dalli']
rhinolopidae = ['Rhinolophus ferrumequinum', 'Rhinolophus sinicus', 'Hipposideros armiger']
moschidae = ['Moschus berezovskii']
indriidae = ['Propithecus coquereli']
leporidae = ['Lepus capensis', 'Bunolagus monticularis', 'Lepus townsendii', 'Lepus timidus', 'Bunolagus monticularis', 'Lepus europaeus', 'Lepus castroviejoi', 'Sylvilagus nutallii', 'Lepus mediterraneus', 'Lepus mandshuricus', 'Lepus coreanus', 'Oryctolagus cuniculus', 'Lepus californicus', 'Lepus saxatilis', 'Pentalagus furnessi']
cheirogaleidae = ['Microcebus murinus']
cercopithecidae = ['Macaca fascicularis', 'Macaca nemestrina', 'Macaca mulatta', 'Macaca fuscata']
dasypodidae = ['Dasypus novemcinctus']
hominidae = ['Pongo pygmaeus', 'Pongo abelii']
sciuridae = ['Sciurus carolinensis']
callitrichidae = ['Leontopithecus chrysomelas']


for seq_record in SeqIO.parse('seqdump.txt', 'fasta'):
    s = ' '
    query = seq_record.id
    query_seq = seq_record.seq
    description = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description:
        description = description.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names = description.split()[1:3]
    names_seq = s.join(names)
     # Verificar se o ID já foi processado
    if names_seq in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq)
        # Adicionando a sequência e a descrição na lista
        lista.append(f'>{names_seq}')
        lista.append(str(query_seq))


    if names_seq in felidae:
        felidae_count += 1
    if names_seq in hyaenidae:
        hyaenidae_count += 1
    if names_seq in delphinidae:
        delphinidae_count += 1
    if names_seq in ursidae:
        ursidae_count += 1
    if names_seq in herpestidae:
        herpestidae_count += 1
    if names_seq in odobenidae:
        odobenidae_count += 1
    if names_seq in phocidae:
        phocidae_count += 1
    if names_seq in otariidae:
        otariidae_count += 1
    if names_seq in canidae:
        canidae_count += 1
    if names_seq in mustelidae:
        mustelidae_count += 1
    if names_seq in physeteridae:
        physeteridae_count += 1
    if names_seq in balaenidae:
        balaenidae_count += 1
    if names_seq in phocoenidae:
        phocoenidae_count += 1
    if names_seq in rhinocerotidae:
        rhinocerotidae_count += 1
    if names_seq in balaenopteridae:
        balaenopteridae_count += 1
    if names_seq in monodontidae:
        monodontidae_count += 1
    if names_seq in pontoporiidae:
        pontoporiidae_count += 1
    if names_seq in lipotidae:
        lipotidae_count += 1
    if names_seq in camelidae:
        camelidae_count += 1
    if names_seq in ziphiidae:
        ziphiidae_count += 1
    if names_seq in hippopotamidae:
        hippopotamidae_count += 1
    if names_seq in vespertilionidae:
        vespertilionidae_count += 1
    if names_seq in suidae:
        suidae_count += 1
    if names_seq in pteropodidae:
        pteropodidae_count += 1
    if names_seq in equidae:
        equidae_count += 1
    if names_seq in cervidae:
        cervidae_count += 1
    if names_seq in tapiridae:
        tapiridae_count += 1
    if names_seq in bovidae:
        bovidae_count += 1
    if names_seq in rhinolopidae:
        rhinolopidae_count += 1
    if names_seq in moschidae:
        moschidae_count += 1
    if names_seq in indriidae:
        indriidae_count += 1
    if names_seq in leporidae:
        leporidae_count += 1
    if names_seq in cheirogaleidae:
        cheirogaleidae_count += 1
    if names_seq in cercopithecidae:
        cercopithecidae_count += 1
    if names_seq in dasypodidae:
        dasypodidae_count += 1
    if names_seq in hominidae:
        hominidae_count += 1
    if names_seq in sciuridae:
        sciuridae_count += 1
    if names_seq in callitrichidae:
        callitrichidae_count += 1
print(f'felidae: {felidae_count}')
print(f'hyaenidae: {hyaenidae_count}')
print(f'delphinidae: {delphinidae_count}')
print(f'ursidae: {ursidae_count}')
print(f'herpestidae: {herpestidae_count}')
print(f'odobenidae:{odobenidae_count}')
print(f'phocidae: {phocidae_count}')
print(f'otariidae: {otariidae_count}')
print(f'canidae: {canidae_count}')
print(f'mustelidae: {mustelidae_count}')
print(f'physeteridae: {physeteridae_count}')
print(f'balaenidae: {balaenidae_count}')
print(f'phocoenidae: {phocoenidae_count}')
print(f'rhinocerotidae: {rhinocerotidae_count}')
print(f'balaenopteridae: {balaenopteridae_count}')
print(f'monodontidae: {monodontidae_count}')
print(f'pontoporiidae: {pontoporiidae_count}')
print(f'lipotidae: {lipotidae_count}')
print(f'camelidae: {camelidae_count}')
print(f'ziphiidae: {ziphiidae_count}')
print(f'hippopotamidae: {hippopotamidae_count}')
print(f'vespertilionidae: {vespertilionidae_count}')
print(f'suidae: {suidae_count}')
print(f'pteropodidae: {pteropodidae_count}')
print(f'equidae: {equidae_count}')
print(f'cervidae: {cervidae_count}')
print(f'tapiridae: {tapiridae_count}')
print(f'bovidae: {bovidae_count}')
print(f'rhinolopidae: {rhinolopidae_count}')
print(f'moschidae: {moschidae_count}')
print(f'indriidae: {indriidae_count}')
print(f'leporidae: {leporidae_count}')
print(f'cheirogaleidae: {cheirogaleidae_count}')
print(f'cercopithecidae: {cercopithecidae_count}')
print(f'dasypodidae: {dasypodidae_count}')
print(f'hominidae: {hominidae_count}')
print(f'sciuridae: {sciuridae_count}')
print(f'callitrichidae: {callitrichidae_count}')


total_count_1 = (
    felidae_count + hyaenidae_count + delphinidae_count + ursidae_count + herpestidae_count +
    odobenidae_count + phocidae_count + otariidae_count + canidae_count + mustelidae_count +
    physeteridae_count + balaenidae_count + phocoenidae_count + rhinocerotidae_count +
    balaenopteridae_count + monodontidae_count + pontoporiidae_count + lipotidae_count +
    camelidae_count + ziphiidae_count + hippopotamidae_count + vespertilionidae_count +
    suidae_count + pteropodidae_count + equidae_count + cervidae_count + tapiridae_count +
    bovidae_count + rhinolopidae_count + moschidae_count + indriidae_count + leporidae_count +
    cheirogaleidae_count + cercopithecidae_count + dasypodidae_count + hominidae_count + sciuridae_count + callitrichidae_count
)

print(f'Total de contagens: {total_count_1}')

# Processamento do arquivo 'seqdump (1).txt'
for seq_record in SeqIO.parse('seqdump (1).txt', 'fasta'):
    s_1 = ' '
    seqdump_1 = seq_record.id
    seqdump_1_seq = seq_record.seq
    description_1 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_1:
        description_1 = description_1.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_1 = description_1.split()[1:3]
    names_seq_1 = s_1.join(names_1)
    
    # Adicionando a sequência e a descrição na lista, se ainda não estiver presente
      # Verificar se o ID já foi processado
    if names_seq_1 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_1)
        lista.append(f'>{names_seq_1}')
        lista.append(str(seqdump_1_seq))

    if names_seq_1 in felidae:
            felidae_count += 1
    if names_seq_1 in hyaenidae:
            hyaenidae_count += 1
    if names_seq_1 in delphinidae:
            delphinidae_count += 1
    if names_seq_1 in ursidae:
            ursidae_count += 1
    if names_seq_1 in herpestidae:
        herpestidae_count += 1
    if names_seq_1 in odobenidae:
        odobenidae_count += 1
    if names_seq_1 in phocidae:
        phocidae_count += 1
    if names_seq_1 in otariidae:
        otariidae_count += 1
    if names_seq_1 in canidae:
        canidae_count += 1
    if names_seq_1 in mustelidae:
        mustelidae_count += 1
    if names_seq_1 in physeteridae:
        physeteridae_count += 1
    if names_seq_1 in balaenidae:
        balaenidae_count += 1
    if names_seq_1 in phocoenidae:
        phocoenidae_count += 1
    if names_seq_1 in rhinocerotidae:
        rhinocerotidae_count += 1
    if names_seq_1 in balaenopteridae:
        balaenopteridae_count += 1
    if names_seq_1 in monodontidae:
        monodontidae_count += 1
    if names_seq_1 in pontoporiidae:
        pontoporiidae_count += 1
    if names_seq_1 in lipotidae:
        lipotidae_count += 1
    if names_seq_1 in camelidae:
        camelidae_count += 1
    if names_seq_1 in ziphiidae:
        ziphiidae_count += 1
    if names_seq_1 in hippopotamidae:
        hippopotamidae_count += 1
    if names_seq_1 in vespertilionidae:
        vespertilionidae_count += 1
    if names_seq_1 in suidae:
        suidae_count += 1
    if names_seq_1 in pteropodidae:
        pteropodidae_count += 1
    if names_seq_1 in equidae:
        equidae_count += 1
    if names_seq_1 in cervidae:
        cervidae_count += 1
    if names_seq_1 in tapiridae:
        tapiridae_count += 1
    if names_seq_1 in bovidae:
        bovidae_count += 1
    if names_seq_1 in rhinolopidae:
        rhinolopidae_count += 1
    if names_seq_1 in moschidae:
        moschidae_count += 1
    if names_seq_1 in indriidae:
        indriidae_count += 1
    if names_seq_1 in leporidae:
        leporidae_count += 1
    if names_seq_1 in cheirogaleidae:
        cheirogaleidae_count += 1
    if names_seq_1 in cercopithecidae:
        cercopithecidae_count += 1
    if names_seq_1 in dasypodidae:
        dasypodidae_count += 1
    if names_seq_1 in hominidae:
        hominidae_count += 1
    if names_seq_1 in sciuridae:
        sciuridae_count += 1
    if names_seq_1 in callitrichidae:
        callitrichidae_count += 1
print(f'felidae: {felidae_count}')
print(f'hyaenidae: {hyaenidae_count}')
print(f'delphinidae: {delphinidae_count}')
print(f'ursidae: {ursidae_count}')
print(f'herpestidae: {herpestidae_count}')
print(f'odobenidae:{odobenidae_count}')
print(f'phocidae: {phocidae_count}')
print(f'otariidae: {otariidae_count}')
print(f'canidae: {canidae_count}')
print(f'mustelidae: {mustelidae_count}')
print(f'physeteridae: {physeteridae_count}')
print(f'balaenidae: {balaenidae_count}')
print(f'phocoenidae: {phocoenidae_count}')
print(f'rhinocerotidae: {rhinocerotidae_count}')
print(f'balaenopteridae: {balaenopteridae_count}')
print(f'monodontidae: {monodontidae_count}')
print(f'pontoporiidae: {pontoporiidae_count}')
print(f'lipotidae: {lipotidae_count}')
print(f'camelidae: {camelidae_count}')
print(f'ziphiidae: {ziphiidae_count}')
print(f'hippopotamidae: {hippopotamidae_count}')
print(f'vespertilionidae: {vespertilionidae_count}')
print(f'suidae: {suidae_count}')
print(f'pteropodidae: {pteropodidae_count}')
print(f'equidae: {equidae_count}')
print(f'cervidae: {cervidae_count}')
print(f'tapiridae: {tapiridae_count}')
print(f'bovidae: {bovidae_count}')
print(f'rhinolopidae: {rhinolopidae_count}')
print(f'moschidae: {moschidae_count}')
print(f'indriidae: {indriidae_count}')
print(f'leporidae: {leporidae_count}')
print(f'cheirogaleidae: {cheirogaleidae_count}')
print(f'cercopithecidae: {cercopithecidae_count}')
print(f'dasypodidae: {dasypodidae_count}')
print(f'hominidae: {hominidae_count}')
print(f'sciuridae: {sciuridae_count}')
print(f'callitrichidae: {callitrichidae_count}')


total_count_2 = (
    felidae_count + hyaenidae_count + delphinidae_count + ursidae_count + herpestidae_count +
    odobenidae_count + phocidae_count + otariidae_count + canidae_count + mustelidae_count +
    physeteridae_count + balaenidae_count + phocoenidae_count + rhinocerotidae_count +
    balaenopteridae_count + monodontidae_count + pontoporiidae_count + lipotidae_count +
    camelidae_count + ziphiidae_count + hippopotamidae_count + vespertilionidae_count +
    suidae_count + pteropodidae_count + equidae_count + cervidae_count + tapiridae_count +
    bovidae_count + rhinolopidae_count + moschidae_count + indriidae_count + leporidae_count +
    cheirogaleidae_count + cercopithecidae_count + dasypodidae_count + hominidae_count + sciuridae_count + callitrichidae_count
)

print(f'Total de contagens: {total_count_2}')



for seq_record in SeqIO.parse('seqdump (2).txt', 'fasta'):
        s_2 = ' '
        seqdump_2 = seq_record.id
        seqdump_2_seq = seq_record.seq
        
        description_2 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
        if 'PREDICTED:' in description_2:
            description_2 = description_2.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
        names_2 = description_2.split()[1:3]
        names_seq_2 = s_2.join(names_2)
          # Verificar se o ID já foi processado
        if names_seq_2 in unique_ids:
            continue  # Ignorar duplicados
        else:
            unique_ids.add(names_seq_2)
            lista.append(f'>{names_seq_2}')
            lista.append(str(seqdump_2_seq))
        if names_seq_2 in felidae:
            felidae_count += 1
        if names_seq_2 in hyaenidae:
            hyaenidae_count += 1
        if names_seq_2 in delphinidae:
            delphinidae_count += 1
        if names_seq_2 in ursidae:
            ursidae_count += 1
        if names_seq_2 in herpestidae:
            herpestidae_count += 1
        if names_seq_2 in odobenidae:
            odobenidae_count += 1
        if names_seq_2 in phocidae:
            phocidae_count += 1
        if names_seq_2 in canidae:
            canidae_count += 1
        if names_seq_2 in mustelidae:
            mustelidae_count += 1
        if names_seq_2 in physeteridae:
            physeteridae_count += 1
        if names_seq_2 in balaenidae:
            balaenidae_count += 1
        if names_seq_2 in phocoenidae:
            phocoenidae_count += 1
        if names_seq_2 in rhinocerotidae:
            rhinocerotidae_count += 1
        if names_seq_2 in balaenopteridae:
            balaenopteridae_count += 1
        if names_seq_2 in monodontidae:
            monodontidae_count += 1
        if names_seq_2 in pontoporiidae:
            pontoporiidae_count += 1
        if names_seq_2 in lipotidae:
            lipotidae_count += 1
        if names_seq_2 in camelidae:
            camelidae_count += 1
        if names_seq_2 in ziphiidae:
            ziphiidae_count += 1
        if names_seq_2 in hippopotamidae:
            hippopotamidae_count += 1
        if names_seq_2 in vespertilionidae:
            vespertilionidae_count += 1
        if names_seq_2 in suidae:
            suidae_count += 1
        if names_seq_2 in pteropodidae:
            pteropodidae_count += 1
        if names_seq_2 in equidae:
            equidae_count += 1
        if names_seq_2 in cervidae:
            cervidae_count += 1
        if names_seq_2 in tapiridae:
            tapiridae_count += 1
        if names_seq_2 in bovidae:
            bovidae_count += 1
        if names_seq_2 in rhinolopidae:
            rhinolopidae_count += 1
        if names_seq_2 in moschidae:
            moschidae_count += 1
        if names_seq_2 in indriidae:
            indriidae_count += 1
        if names_seq_2 in leporidae:
            leporidae_count += 1
        if names_seq_2 in cheirogaleidae:
            cheirogaleidae_count += 1
        if names_seq_2 in cercopithecidae:
            cercopithecidae_count += 1
        if names_seq_2 in dasypodidae:
            dasypodidae_count += 1
        if names_seq_2 in hominidae:
            hominidae_count += 1
        if names_seq_2 in sciuridae:
            sciuridae_count += 1
        if names_seq_2 in callitrichidae:
            callitrichidae_count += 1
print(f'felidae: {felidae_count}')
print(f'hyaenidae: {hyaenidae_count}')
print(f'delphinidae: {delphinidae_count}')
print(f'ursidae: {ursidae_count}')
print(f'herpestidae: {herpestidae_count}')
print(f'odobenidae:{odobenidae_count}')
print(f'phocidae: {phocidae_count}')
print(f'otariidae: {otariidae_count}')
print(f'canidae: {canidae_count}')
print(f'mustelidae: {mustelidae_count}')
print(f'physeteridae: {physeteridae_count}')
print(f'balaenidae: {balaenidae_count}')
print(f'phocoenidae: {phocoenidae_count}')
print(f'rhinocerotidae: {rhinocerotidae_count}')
print(f'balaenopteridae: {balaenopteridae_count}')
print(f'monodontidae: {monodontidae_count}')
print(f'pontoporiidae: {pontoporiidae_count}')
print(f'lipotidae: {lipotidae_count}')
print(f'camelidae: {camelidae_count}')
print(f'ziphiidae: {ziphiidae_count}')
print(f'hippopotamidae: {hippopotamidae_count}')
print(f'vespertilionidae: {vespertilionidae_count}')
print(f'suidae: {suidae_count}')
print(f'pteropodidae: {pteropodidae_count}')
print(f'equidae: {equidae_count}')
print(f'cervidae: {cervidae_count}')
print(f'tapiridae: {tapiridae_count}')
print(f'bovidae: {bovidae_count}')
print(f'rhinolopidae: {rhinolopidae_count}')
print(f'moschidae: {moschidae_count}')
print(f'indriidae: {indriidae_count}')
print(f'leporidae: {leporidae_count}')
print(f'cheirogaleidae: {cheirogaleidae_count}')
print(f'cercopithecidae: {cercopithecidae_count}')
print(f'dasypodidae: {dasypodidae_count}')
print(f'hominidae: {hominidae_count}')
print(f'sciuridae: {sciuridae_count}')
print(f'callitrichidae: {callitrichidae_count}')


total_count_7 = (
    felidae_count + hyaenidae_count + delphinidae_count + ursidae_count + herpestidae_count +
    odobenidae_count + phocidae_count + otariidae_count + canidae_count + mustelidae_count +
    physeteridae_count + balaenidae_count + phocoenidae_count + rhinocerotidae_count +
    balaenopteridae_count + monodontidae_count + pontoporiidae_count + lipotidae_count +
    camelidae_count + ziphiidae_count + hippopotamidae_count + vespertilionidae_count +
    suidae_count + pteropodidae_count + equidae_count + cervidae_count + tapiridae_count +
    bovidae_count + rhinolopidae_count + moschidae_count + indriidae_count + leporidae_count +
    cheirogaleidae_count + cercopithecidae_count + dasypodidae_count + hominidae_count + sciuridae_count + callitrichidae_count
)

print(f'Total de contagens: {total_count_7}')

for seq_record in SeqIO.parse('seqdump (3).txt', 'fasta'):
        seqdump_3 = seq_record.id
        seqdump_3_seq = seq_record.seq
        s_3 = ' '
        description_3 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
        if 'PREDICTED:' in description_3:
            description_3 = description_3.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
        names_3 = description_3.split()[1:3]
        names_seq_3 = s_3.join(names_3)
        if names_seq_3 in unique_ids:
            continue  # Ignorar duplicados
        else:
            unique_ids.add(names_seq_3)
            lista.append(f'>{names_seq_3}')
            lista.append(str(seqdump_3_seq))
for seq_record in SeqIO.parse('seqdump (4).txt', 'fasta'):
    s_4 = ' '
    seqdump_4 = seq_record.id
    seqdump_4_seq = seq_record.seq
    description_4 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_4:
        description_4 = description_4.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_4 = description_4.split()[1:3]
    names_seq_4 = s_4.join(names_4)
    
    if names_seq_4 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_4)
        lista.append(f'>{names_seq_4}')
        lista.append(str(seqdump_4_seq))

for seq_record in SeqIO.parse('seqdump (5).txt', 'fasta'):
    s_5 = ' '
    seqdump_5 = seq_record.id
    seqdump_5_seq = seq_record.seq
    description_5 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_5:
        description_5 = description_5.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_5 = description_5.split()[1:3]
    names_seq_5 = s_5.join(names_5)
    
    # Adicionando a sequência e a descrição na lista, se ainda não estiver presente
    if names_seq_5 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_5)
        lista.append(f'>{names_seq_5}')
        lista.append(str(seqdump_5_seq))
for seq_record in SeqIO.parse('seqdump (6).txt', 'fasta'):
    s_6 = ' '
    seqdump_6 = seq_record.id
    seqdump_6_seq = seq_record.seq
    description_6 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_6:
        description_6 = description_6.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_6 = description_6.split()[1:3]
    names_seq_6 = s_6.join(names_6)
    
    if names_seq_6 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_6)
        lista.append(f'>{names_seq_6}')
        lista.append(str(seqdump_6_seq))
for seq_record in SeqIO.parse('seqdump (7).txt', 'fasta'):
    s_7 = ' '
    seqdump_7 = seq_record.id
    seqdump_7_seq = seq_record.seq
    description_7 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_7:
        description_7 = description_7.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_7 = description_7.split()[1:3]
    names_seq_7 = s_7.join(names_7)
    
    # Adicionando a sequência e a descrição na lista, se ainda não estiver presente
    if names_seq_7 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_7)
        lista.append(f'>{names_seq_7}')
        lista.append(str(seqdump_7_seq))
for seq_record in SeqIO.parse('seqdump (8).txt', 'fasta'):
    s_8 = ' '
    seqdump_8 = seq_record.id
    seqdump_8_seq = seq_record.seq
    description_8 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_8:
        description_8 = description_8.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_8 = description_8.split()[1:3]
    names_seq_8 = s_8.join(names_8)
    
    # Adicionando a sequência e a descrição na lista, se ainda não estiver presente
    if names_seq_8 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_8)
        lista.append(f'>{names_seq_8}') # MUDANÇAAAAAAAAAAAAAAAAAAAAAAAAAA
        lista.append(str(seqdump_8_seq))

    if names_seq_8 in felidae:
            felidae_count += 1
    if names_seq_8 in hyaenidae:
            hyaenidae_count += 1
    if names_seq_8 in delphinidae:
            delphinidae_count += 1
    if names_seq_8 in ursidae:
            ursidae_count += 1
    if names_seq_8 in herpestidae:
        herpestidae_count += 1
    if names_seq_8 in odobenidae:
        odobenidae_count += 1
    if names_seq_8 in phocidae:
        phocidae_count += 1
    if names_seq_8 in otariidae:
        otariidae_count += 1
    if names_seq_8 in canidae:
        canidae_count += 1
    if names_seq_8 in mustelidae:
        mustelidae_count += 1
    if names_seq_8 in physeteridae:
        physeteridae_count += 1
    if names_seq_8 in balaenidae:
        balaenidae_count += 1
    if names_seq_8 in phocoenidae:
        phocoenidae_count += 1
    if names_seq_8 in rhinocerotidae:
        rhinocerotidae_count += 1
    if names_seq_8 in balaenopteridae:
        balaenopteridae_count += 1
    if names_seq_8 in monodontidae:
        monodontidae_count += 1
    if names_seq_8 in pontoporiidae:
        pontoporiidae_count += 1
    if names_seq_8 in lipotidae:
        lipotidae_count += 1
    if names_seq_8 in camelidae:
        camelidae_count += 1
    if names_seq_8 in ziphiidae:
        ziphiidae_count += 1
    if names_seq_8 in hippopotamidae:
        hippopotamidae_count += 1
    if names_seq_8 in vespertilionidae:
        vespertilionidae_count += 1
    if names_seq_8 in suidae:
        suidae_count += 1
    if names_seq_8 in pteropodidae:
        pteropodidae_count += 1
    if names_seq_8 in equidae:
        equidae_count += 1
    if names_seq_8 in cervidae:
        cervidae_count += 1
    if names_seq_8 in tapiridae:
        tapiridae_count += 1
    if names_seq_8 in bovidae:
        bovidae_count += 1
    if names_seq_8 in rhinolopidae:
        rhinolopidae_count += 1
    if names_seq_8 in moschidae:
        moschidae_count += 1
    if names_seq_8 in indriidae:
        indriidae_count += 1
    if names_seq_8 in leporidae:
        leporidae_count += 1
    if names_seq_8 in cheirogaleidae:
        cheirogaleidae_count += 1
    if names_seq_8 in cercopithecidae:
        cercopithecidae_count += 1
    if names_seq_8 in dasypodidae:
        dasypodidae_count += 1
    if names_seq_8 in hominidae:
        hominidae_count += 1
    if names_seq_8 in sciuridae:
        sciuridae_count += 1
    if names_seq_8 in callitrichidae:
        callitrichidae_count += 1
print(f'felidae: {felidae_count}')
print(f'hyaenidae: {hyaenidae_count}')
print(f'delphinidae: {delphinidae_count}')
print(f'ursidae: {ursidae_count}')
print(f'herpestidae: {herpestidae_count}')
print(f'odobenidae:{odobenidae_count}')
print(f'phocidae: {phocidae_count}')
print(f'otariidae: {otariidae_count}')
print(f'canidae: {canidae_count}')
print(f'mustelidae: {mustelidae_count}')
print(f'physeteridae: {physeteridae_count}')
print(f'balaenidae: {balaenidae_count}')
print(f'phocoenidae: {phocoenidae_count}')
print(f'rhinocerotidae: {rhinocerotidae_count}')
print(f'balaenopteridae: {balaenopteridae_count}')
print(f'monodontidae: {monodontidae_count}')
print(f'pontoporiidae: {pontoporiidae_count}')
print(f'lipotidae: {lipotidae_count}')
print(f'camelidae: {camelidae_count}')
print(f'ziphiidae: {ziphiidae_count}')
print(f'hippopotamidae: {hippopotamidae_count}')
print(f'vespertilionidae: {vespertilionidae_count}')
print(f'suidae: {suidae_count}')
print(f'pteropodidae: {pteropodidae_count}')
print(f'equidae: {equidae_count}')
print(f'cervidae: {cervidae_count}')
print(f'tapiridae: {tapiridae_count}')
print(f'bovidae: {bovidae_count}')
print(f'rhinolopidae: {rhinolopidae_count}')
print(f'moschidae: {moschidae_count}')
print(f'indriidae: {indriidae_count}')
print(f'leporidae: {leporidae_count}')
print(f'cheirogaleidae: {cheirogaleidae_count}')
print(f'cercopithecidae: {cercopithecidae_count}')
print(f'dasypodidae: {dasypodidae_count}')
print(f'hominidae: {hominidae_count}')
print(f'sciuridae: {sciuridae_count}')
print(f'callitrichidae: {callitrichidae_count}')


total_count_3 = (
    felidae_count + hyaenidae_count + delphinidae_count + ursidae_count + herpestidae_count +
    odobenidae_count + phocidae_count + otariidae_count + canidae_count + mustelidae_count +
    physeteridae_count + balaenidae_count + phocoenidae_count + rhinocerotidae_count +
    balaenopteridae_count + monodontidae_count + pontoporiidae_count + lipotidae_count +
    camelidae_count + ziphiidae_count + hippopotamidae_count + vespertilionidae_count +
    suidae_count + pteropodidae_count + equidae_count + cervidae_count + tapiridae_count +
    bovidae_count + rhinolopidae_count + moschidae_count + indriidae_count + leporidae_count +
    cheirogaleidae_count + cercopithecidae_count + dasypodidae_count + hominidae_count + sciuridae_count + callitrichidae_count
)

print(f'Total de contagens: {total_count_3}')



for seq_record in SeqIO.parse('seqdump (9).txt', 'fasta'):
    s_9 = ' '
    seqdump_9 = seq_record.id
    seqdump_9_seq = seq_record.seq
    description_9 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_9:
        description_9 = description_9.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_9 = description_9.split()[1:3]
    names_seq_9 = s_9.join(names_9)
    
    if names_seq_9 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_9)
        lista.append(f'>{names_seq_9}')
        lista.append(str(seqdump_9_seq))
for seq_record in SeqIO.parse('seqdump (10).txt', 'fasta'):
    s_10 = ' '
    seqdump_10 = seq_record.id
    seqdump_10_seq = seq_record.seq
    description_10 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_10:
        description_10 = description_10.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_10 = description_10.split()[1:3]
    names_seq_10 = s_10.join(names_10)
    
    # Adicionando a sequência e a descrição na lista, se ainda não estiver presente
    if names_seq_10 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_10)
        lista.append(f'>{names_seq_10}')
        lista.append(str(seqdump_10_seq))
for seq_record in SeqIO.parse('seqdump (11).txt', 'fasta'):
    s_11 = ' '
    seqdump_11 = seq_record.id
    seqdump_11_seq = seq_record.seq
    description_11 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_11:
        description_11 = description_11.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_11 = description_11.split()[1:3]
    names_seq_11 = s_11.join(names_11)
    
    # Adicionando a sequência e a descrição na lista, se ainda não estiver presente
    if names_seq_11 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_11)
        lista.append(f'>{names_seq_11}')
        lista.append(str(seqdump_11_seq)) # HOUVE MUDANÇAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA A PARTIR DAQUI
    if names_seq_11 in felidae:
            felidae_count += 1
    if names_seq_11 in hyaenidae:
            hyaenidae_count += 1
    if names_seq_11 in delphinidae:
            delphinidae_count += 1
    if names_seq_11 in ursidae:
            ursidae_count += 1
    if names_seq_11 in herpestidae:
        herpestidae_count += 1
    if names_seq_11 in odobenidae:
        odobenidae_count += 1
    if names_seq_11 in phocidae:
        phocidae_count += 1
    if names_seq_11 in otariidae:
        otariidae_count += 1
    if names_seq_11 in canidae:
        canidae_count += 1
    if names_seq_11 in mustelidae:
        mustelidae_count += 1
    if names_seq_11 in physeteridae:
        physeteridae_count += 1
    if names_seq_11 in balaenidae:
        balaenidae_count += 1
    if names_seq_11 in phocoenidae:
        phocoenidae_count += 1
    if names_seq_11 in rhinocerotidae:
        rhinocerotidae_count += 1
    if names_seq_11 in balaenopteridae:
        balaenopteridae_count += 1
    if names_seq_11 in monodontidae:
        monodontidae_count += 1
    if names_seq_11 in pontoporiidae:
        pontoporiidae_count += 1
    if names_seq_11 in lipotidae:
        lipotidae_count += 1
    if names_seq_11 in camelidae:
        camelidae_count += 1
    if names_seq_11 in ziphiidae:
        ziphiidae_count += 1
    if names_seq_11 in hippopotamidae:
        hippopotamidae_count += 1
    if names_seq_11 in vespertilionidae:
        vespertilionidae_count += 1
    if names_seq_11 in suidae:
        suidae_count += 1
    if names_seq_11 in pteropodidae:
        pteropodidae_count += 1
    if names_seq_11 in equidae:
        equidae_count += 1
    if names_seq_11 in cervidae:
        cervidae_count += 1
    if names_seq_11 in tapiridae:
        tapiridae_count += 1
    if names_seq_11 in bovidae:
        bovidae_count += 1
    if names_seq_11 in rhinolopidae:
        rhinolopidae_count += 1
    if names_seq_11 in moschidae:
        moschidae_count += 1
    if names_seq_11 in indriidae:
        indriidae_count += 1
    if names_seq_11 in leporidae:
        leporidae_count += 1
    if names_seq_11 in cheirogaleidae:
        cheirogaleidae_count += 1
    if names_seq_11 in cercopithecidae:
        cercopithecidae_count += 1
    if names_seq_11 in dasypodidae:
        dasypodidae_count += 1
    if names_seq_11 in hominidae:
        hominidae_count += 1
    if names_seq_11 in sciuridae:
        sciuridae_count += 1
    if names_seq_11 in callitrichidae:
        callitrichidae_count += 1
print(f'felidae: {felidae_count}')
print(f'hyaenidae: {hyaenidae_count}')
print(f'delphinidae: {delphinidae_count}')
print(f'ursidae: {ursidae_count}')
print(f'herpestidae: {herpestidae_count}')
print(f'odobenidae:{odobenidae_count}')
print(f'phocidae: {phocidae_count}')
print(f'otariidae: {otariidae_count}')
print(f'canidae: {canidae_count}')
print(f'mustelidae: {mustelidae_count}')
print(f'physeteridae: {physeteridae_count}')
print(f'balaenidae: {balaenidae_count}')
print(f'phocoenidae: {phocoenidae_count}')
print(f'rhinocerotidae: {rhinocerotidae_count}')
print(f'balaenopteridae: {balaenopteridae_count}')
print(f'monodontidae: {monodontidae_count}')
print(f'pontoporiidae: {pontoporiidae_count}')
print(f'lipotidae: {lipotidae_count}')
print(f'camelidae: {camelidae_count}')
print(f'ziphiidae: {ziphiidae_count}')
print(f'hippopotamidae: {hippopotamidae_count}')
print(f'vespertilionidae: {vespertilionidae_count}')
print(f'suidae: {suidae_count}')
print(f'pteropodidae: {pteropodidae_count}')
print(f'equidae: {equidae_count}')
print(f'cervidae: {cervidae_count}')
print(f'tapiridae: {tapiridae_count}')
print(f'bovidae: {bovidae_count}')
print(f'rhinolopidae: {rhinolopidae_count}')
print(f'moschidae: {moschidae_count}')
print(f'indriidae: {indriidae_count}')
print(f'leporidae: {leporidae_count}')
print(f'cheirogaleidae: {cheirogaleidae_count}')
print(f'cercopithecidae: {cercopithecidae_count}')
print(f'dasypodidae: {dasypodidae_count}')
print(f'hominidae: {hominidae_count}')
print(f'sciuridae: {sciuridae_count}')
print(f'callitrichidae: {callitrichidae_count}')


total_count_3 = (
    felidae_count + hyaenidae_count + delphinidae_count + ursidae_count + herpestidae_count +
    odobenidae_count + phocidae_count + otariidae_count + canidae_count + mustelidae_count +
    physeteridae_count + balaenidae_count + phocoenidae_count + rhinocerotidae_count +
    balaenopteridae_count + monodontidae_count + pontoporiidae_count + lipotidae_count +
    camelidae_count + ziphiidae_count + hippopotamidae_count + vespertilionidae_count +
    suidae_count + pteropodidae_count + equidae_count + cervidae_count + tapiridae_count +
    bovidae_count + rhinolopidae_count + moschidae_count + indriidae_count + leporidae_count +
    cheirogaleidae_count + cercopithecidae_count + dasypodidae_count + hominidae_count + sciuridae_count + callitrichidae_count
)

print(f'Total de contagens: {total_count_3}')
        
for seq_record in SeqIO.parse('seqdump (12).txt', 'fasta'):
    s_12 = ' '
    seqdump_12 = seq_record.id
    seqdump_12_seq = seq_record.seq
    description_12 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_12:
        description_12 = description_12.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_12 = description_12.split()[1:3]
    names_seq_12 = s_12.join(names_12)
    
    # Adicionando a sequência e a descrição na lista, se ainda não estiver presente
    if names_seq_12 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_12)
        lista.append(f'>{names_seq_12}')
        lista.append(str(seqdump_12_seq))
for seq_record in SeqIO.parse('seqdump (13).txt', 'fasta'):
    s_13 = ' '
    seqdump_13 = seq_record.id
    seqdump_13_seq = seq_record.seq
    description_13 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_13:
        description_13 = description_13.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_13 = description_13.split()[1:3]
    names_seq_13 = s_13.join(names_13)
    
    if names_seq_13 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_13)
        lista.append(f'>{names_seq_13}')
        lista.append(str(seqdump_13_seq))
for seq_record in SeqIO.parse('seqdump (14).txt', 'fasta'):
    s_14 = ' '
    seqdump_14 = seq_record.id
    seqdump_14_seq = seq_record.seq
    description_14 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_14:
        description_14 = description_14.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_14 = description_14.split()[1:3]
    names_seq_14 = s_14.join(names_14)
    
    if names_seq_14 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_14)
        lista.append(f'>{names_seq_14}')
        lista.append(str(seqdump_14_seq))
for seq_record in SeqIO.parse('seqdump (15).txt', 'fasta'):
    s_15 = ' '
    seqdump_15 = seq_record.id
    seqdump_15_seq = seq_record.seq
    description_15 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_15:
        description_15 = description_15.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_15 = description_15.split()[1:3]
    names_seq_15 = s_15.join(names_15)
    
    if names_seq_15 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_15)
        lista.append(f'{names_seq_15}')
for seq_record in SeqIO.parse('seqdump (16).txt', 'fasta'):
    s_16 = ' '
    seqdump_16 = seq_record.id
    seqdump_16_seq = seq_record.seq
    description_16 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_16:
        description_16 = description_16.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_16 = description_16.split()[1:3]
    names_seq_16 = s_16.join(names_16)
    
    if names_seq_16 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_16)
        lista.append(f'>{names_seq_16}')
        lista.append(str(seqdump_16_seq))
for seq_record in SeqIO.parse('seqdump (17).txt', 'fasta'):
    s_17 = ' '
    seqdump_17 = seq_record.id
    seqdump_17_seq = seq_record.seq
    description_17 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_17:
        description_17 = description_17.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_17 = description_17.split()[1:3]
    names_seq_17 = s_17.join(names_17)
    
    if names_seq_17 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_17)
        lista.append(f'>{names_seq_17}') #OVIS DALLI MAIS UMM AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        lista.append(str(seqdump_17_seq))
    if names_seq_17 in felidae:
            felidae_count += 1
    if names_seq_17 in hyaenidae:
            hyaenidae_count += 1
    if names_seq_17 in delphinidae:
            delphinidae_count += 1
    if names_seq_17 in ursidae:
            ursidae_count += 1
    if names_seq_17 in herpestidae:
        herpestidae_count += 1
    if names_seq_17 in odobenidae:
        odobenidae_count += 1
    if names_seq_17 in phocidae:
        phocidae_count += 1
    if names_seq_17 in otariidae:
        otariidae_count += 1
    if names_seq_17 in canidae:
        canidae_count += 1
    if names_seq_17 in mustelidae:
        mustelidae_count += 1
    if names_seq_17 in physeteridae:
        physeteridae_count += 1
    if names_seq_17 in balaenidae:
        balaenidae_count += 1
    if names_seq_17 in phocoenidae:
        phocoenidae_count += 1
    if names_seq_17 in rhinocerotidae:
        rhinocerotidae_count += 1
    if names_seq_17 in balaenopteridae:
        balaenopteridae_count += 1
    if names_seq_17 in monodontidae:
        monodontidae_count += 1
    if names_seq_17 in pontoporiidae:
        pontoporiidae_count += 1
    if names_seq_17 in lipotidae:
        lipotidae_count += 1
    if names_seq_17 in camelidae:
        camelidae_count += 1
    if names_seq_17 in ziphiidae:
        ziphiidae_count += 1
    if names_seq_17 in hippopotamidae:
        hippopotamidae_count += 1
    if names_seq_17 in vespertilionidae:
        vespertilionidae_count += 1
    if names_seq_17 in suidae:
        suidae_count += 1
    if names_seq_17 in pteropodidae:
        pteropodidae_count += 1
    if names_seq_17 in equidae:
        equidae_count += 1
    if names_seq_17 in cervidae:
        cervidae_count += 1
    if names_seq_17 in tapiridae:
        tapiridae_count += 1
    if names_seq_17 in bovidae:
        bovidae_count += 1
    if names_seq_17 in rhinolopidae:
        rhinolopidae_count += 1
    if names_seq_17 in moschidae:
        moschidae_count += 1
    if names_seq_17 in indriidae:
        indriidae_count += 1
    if names_seq_17 in leporidae:
        leporidae_count += 1
    if names_seq_17 in cheirogaleidae:
        cheirogaleidae_count += 1
    if names_seq_17 in cercopithecidae:
        cercopithecidae_count += 1
    if names_seq_17 in dasypodidae:
        dasypodidae_count += 1
    if names_seq_17 in hominidae:
        hominidae_count += 1
    if names_seq_17 in sciuridae:
        sciuridae_count += 1
    if names_seq_17 in callitrichidae:
        callitrichidae_count += 1
print(f'felidae: {felidae_count}')
print(f'hyaenidae: {hyaenidae_count}')
print(f'delphinidae: {delphinidae_count}')
print(f'ursidae: {ursidae_count}')
print(f'herpestidae: {herpestidae_count}')
print(f'odobenidae:{odobenidae_count}')
print(f'phocidae: {phocidae_count}')
print(f'otariidae: {otariidae_count}')
print(f'canidae: {canidae_count}')
print(f'mustelidae: {mustelidae_count}')
print(f'physeteridae: {physeteridae_count}')
print(f'balaenidae: {balaenidae_count}')
print(f'phocoenidae: {phocoenidae_count}')
print(f'rhinocerotidae: {rhinocerotidae_count}')
print(f'balaenopteridae: {balaenopteridae_count}')
print(f'monodontidae: {monodontidae_count}')
print(f'pontoporiidae: {pontoporiidae_count}')
print(f'lipotidae: {lipotidae_count}')
print(f'camelidae: {camelidae_count}')
print(f'ziphiidae: {ziphiidae_count}')
print(f'hippopotamidae: {hippopotamidae_count}')
print(f'vespertilionidae: {vespertilionidae_count}')
print(f'suidae: {suidae_count}')
print(f'pteropodidae: {pteropodidae_count}')
print(f'equidae: {equidae_count}')
print(f'cervidae: {cervidae_count}')
print(f'tapiridae: {tapiridae_count}')
print(f'bovidae: {bovidae_count}')
print(f'rhinolopidae: {rhinolopidae_count}')
print(f'moschidae: {moschidae_count}')
print(f'indriidae: {indriidae_count}')
print(f'leporidae: {leporidae_count}')
print(f'cheirogaleidae: {cheirogaleidae_count}')
print(f'cercopithecidae: {cercopithecidae_count}')
print(f'dasypodidae: {dasypodidae_count}')
print(f'hominidae: {hominidae_count}')
print(f'sciuridae: {sciuridae_count}')
print(f'callitrichidae: {callitrichidae_count}')


total_count_4 = (
    felidae_count + hyaenidae_count + delphinidae_count + ursidae_count + herpestidae_count +
    odobenidae_count + phocidae_count + otariidae_count + canidae_count + mustelidae_count +
    physeteridae_count + balaenidae_count + phocoenidae_count + rhinocerotidae_count +
    balaenopteridae_count + monodontidae_count + pontoporiidae_count + lipotidae_count +
    camelidae_count + ziphiidae_count + hippopotamidae_count + vespertilionidae_count +
    suidae_count + pteropodidae_count + equidae_count + cervidae_count + tapiridae_count +
    bovidae_count + rhinolopidae_count + moschidae_count + indriidae_count + leporidae_count +
    cheirogaleidae_count + cercopithecidae_count + dasypodidae_count + hominidae_count + sciuridae_count + callitrichidae_count
)

print(f'Total de contagens: {total_count_4}')

for seq_record in SeqIO.parse('seqdump (18).txt', 'fasta'):
    s_18 = ' '
    seqdump_18 = seq_record.id
    seqdump_18_seq = seq_record.seq
    description_18 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_18:
        description_18 = description_18.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_18 = description_18.split()[1:3]
    names_seq_18 = s.join(names_18)
    
    if names_seq_18 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_18)
        lista.append(f'>{names_seq_18}') # Lepus saxatilis e Lepus californicus AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        lista.append(str(seqdump_18_seq))

    if names_seq_18 in felidae:
        felidae_count += 1
    if names_seq_18 in hyaenidae:
        hyaenidae_count += 1
    if names_seq_18 in delphinidae:
        delphinidae_count += 1
    if names_seq_18 in ursidae:
        ursidae_count += 1
    if names_seq_18 in herpestidae:
        herpestidae_count += 1
    if names_seq_18 in odobenidae:
        odobenidae_count += 1
    if names_seq_18 in phocidae:
        phocidae_count += 1
    if names_seq_18 in otariidae:
        otariidae_count += 1
    if names_seq_18 in canidae:
        canidae_count += 1
    if names_seq_18 in mustelidae:
        mustelidae_count += 1
    if names_seq_18 in physeteridae:
        physeteridae_count += 1
    if names_seq_18 in balaenidae:
        balaenidae_count += 1
    if names_seq_18 in phocoenidae:
        phocoenidae_count += 1
    if names_seq_18 in rhinocerotidae:
        rhinocerotidae_count += 1
    if names_seq_18 in balaenopteridae:
        balaenopteridae_count += 1
    if names_seq_18 in monodontidae:
        monodontidae_count += 1
    if names_seq_18 in pontoporiidae:
        pontoporiidae_count += 1
    if names_seq_18 in lipotidae:
        lipotidae_count += 1
    if names_seq_18 in camelidae:
        camelidae_count += 1
    if names_seq_18 in ziphiidae:
        ziphiidae_count += 1
    if names_seq_18 in hippopotamidae:
        hippopotamidae_count += 1
    if names_seq_18 in vespertilionidae:
        vespertilionidae_count += 1
    if names_seq_18 in suidae:
        suidae_count += 1
    if names_seq_18 in pteropodidae:
        pteropodidae_count += 1
    if names_seq_18 in equidae:
        equidae_count += 1
    if names_seq_18 in cervidae:
        cervidae_count += 1
    if names_seq_18 in tapiridae:
        tapiridae_count += 1
    if names_seq_18 in bovidae:
        bovidae_count += 1
    if names_seq_18 in rhinolopidae:
        rhinolopidae_count += 1
    if names_seq_18 in moschidae:
        moschidae_count += 1
    if names_seq_18 in indriidae:
        indriidae_count += 1
    if names_seq_18 in leporidae:
        leporidae_count += 1
    if names_seq_18 in cheirogaleidae:
        cheirogaleidae_count += 1
    if names_seq_18 in cercopithecidae:
        cercopithecidae_count += 1
    if names_seq_18 in dasypodidae:
        dasypodidae_count += 1
    if names_seq_18 in hominidae:
        hominidae_count += 1
    if names_seq_18 in sciuridae:
        sciuridae_count += 1
    if names_seq_18 in callitrichidae:
        callitrichidae_count += 1
print(f'felidae: {felidae_count}')
print(f'hyaenidae: {hyaenidae_count}')
print(f'delphinidae: {delphinidae_count}')
print(f'ursidae: {ursidae_count}')
print(f'herpestidae: {herpestidae_count}')
print(f'odobenidae:{odobenidae_count}')
print(f'phocidae: {phocidae_count}')
print(f'otariidae: {otariidae_count}')
print(f'canidae: {canidae_count}')
print(f'mustelidae: {mustelidae_count}')
print(f'physeteridae: {physeteridae_count}')
print(f'balaenidae: {balaenidae_count}')
print(f'phocoenidae: {phocoenidae_count}')
print(f'rhinocerotidae: {rhinocerotidae_count}')
print(f'balaenopteridae: {balaenopteridae_count}')
print(f'monodontidae: {monodontidae_count}')
print(f'pontoporiidae: {pontoporiidae_count}')
print(f'lipotidae: {lipotidae_count}')
print(f'camelidae: {camelidae_count}')
print(f'ziphiidae: {ziphiidae_count}')
print(f'hippopotamidae: {hippopotamidae_count}')
print(f'vespertilionidae: {vespertilionidae_count}')
print(f'suidae: {suidae_count}')
print(f'pteropodidae: {pteropodidae_count}')
print(f'equidae: {equidae_count}')
print(f'cervidae: {cervidae_count}')
print(f'tapiridae: {tapiridae_count}')
print(f'bovidae: {bovidae_count}')
print(f'rhinolopidae: {rhinolopidae_count}')
print(f'moschidae: {moschidae_count}')
print(f'indriidae: {indriidae_count}')
print(f'leporidae: {leporidae_count}')
print(f'cheirogaleidae: {cheirogaleidae_count}')
print(f'cercopithecidae: {cercopithecidae_count}')
print(f'dasypodidae: {dasypodidae_count}')
print(f'hominidae: {hominidae_count}')
print(f'sciuridae: {sciuridae_count}')
print(f'callitrichidae: {callitrichidae_count}')


total_count_8 = (
    felidae_count + hyaenidae_count + delphinidae_count + ursidae_count + herpestidae_count +
    odobenidae_count + phocidae_count + otariidae_count + canidae_count + mustelidae_count +
    physeteridae_count + balaenidae_count + phocoenidae_count + rhinocerotidae_count +
    balaenopteridae_count + monodontidae_count + pontoporiidae_count + lipotidae_count +
    camelidae_count + ziphiidae_count + hippopotamidae_count + vespertilionidae_count +
    suidae_count + pteropodidae_count + equidae_count + cervidae_count + tapiridae_count +
    bovidae_count + rhinolopidae_count + moschidae_count + indriidae_count + leporidae_count +
    cheirogaleidae_count + cercopithecidae_count + dasypodidae_count + hominidae_count + sciuridae_count + callitrichidae_count
)

print(f'Total de contagens: {total_count_8}')        
for seq_record in SeqIO.parse('seqdump (19).txt', 'fasta'):
    s_19 = ' '
    seqdump_19 = seq_record.id
    seqdump_19_seq = seq_record.seq
    description_19 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_19:
        description_19 = description_19.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_19 = description_19.split()[1:3]
    names_seq_19 = s_19.join(names_19)
    
    if names_seq_19 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_19)
        lista.append(f'>{names_seq_19}') # Pentalagus furnessi MAIS UM AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
        lista.append(str(seqdump_19_seq))
    if names_seq_19 in felidae:
            felidae_count += 1
    if names_seq_19 in hyaenidae:
            hyaenidae_count += 1
    if names_seq_19 in delphinidae:
            delphinidae_count += 1
    if names_seq_19 in ursidae:
            ursidae_count += 1
    if names_seq_19 in herpestidae:
        herpestidae_count += 1
    if names_seq_19 in odobenidae:
        odobenidae_count += 1
    if names_seq_19 in phocidae:
        phocidae_count += 1
    if names_seq_1 in otariidae:
        otariidae_count += 1
    if names_seq_19 in canidae:
        canidae_count += 1
    if names_seq_19 in mustelidae:
        mustelidae_count += 1
    if names_seq_19 in physeteridae:
        physeteridae_count += 1
    if names_seq_19 in balaenidae:
        balaenidae_count += 1
    if names_seq_19 in phocoenidae:
        phocoenidae_count += 1
    if names_seq_19 in rhinocerotidae:
        rhinocerotidae_count += 1
    if names_seq_19 in balaenopteridae:
        balaenopteridae_count += 1
    if names_seq_19 in monodontidae:
        monodontidae_count += 1
    if names_seq_19 in pontoporiidae:
        pontoporiidae_count += 1
    if names_seq_19 in lipotidae:
        lipotidae_count += 1
    if names_seq_19 in camelidae:
        camelidae_count += 1
    if names_seq_19 in ziphiidae:
        ziphiidae_count += 1
    if names_seq_19 in hippopotamidae:
        hippopotamidae_count += 1
    if names_seq_19 in vespertilionidae:
        vespertilionidae_count += 1
    if names_seq_19 in suidae:
        suidae_count += 1
    if names_seq_19 in pteropodidae:
        pteropodidae_count += 1
    if names_seq_19 in equidae:
        equidae_count += 1
    if names_seq_19 in cervidae:
        cervidae_count += 1
    if names_seq_19 in tapiridae:
        tapiridae_count += 1
    if names_seq_19 in bovidae:
        bovidae_count += 1
    if names_seq_19 in rhinolopidae:
        rhinolopidae_count += 1
    if names_seq_19 in moschidae:
        moschidae_count += 1
    if names_seq_19 in indriidae:
        indriidae_count += 1
    if names_seq_19 in leporidae:
        leporidae_count += 1
    if names_seq_19 in cheirogaleidae:
        cheirogaleidae_count += 1
    if names_seq_19 in cercopithecidae:
        cercopithecidae_count += 1
    if names_seq_19 in dasypodidae:
        dasypodidae_count += 1
    if names_seq_19 in hominidae:
        hominidae_count += 1
    if names_seq_19 in sciuridae:
        sciuridae_count += 1
    if names_seq_19 in callitrichidae:
        callitrichidae_count += 1
print(f'felidae: {felidae_count}')
print(f'hyaenidae: {hyaenidae_count}')
print(f'delphinidae: {delphinidae_count}')
print(f'ursidae: {ursidae_count}')
print(f'herpestidae: {herpestidae_count}')
print(f'odobenidae:{odobenidae_count}')
print(f'phocidae: {phocidae_count}')
print(f'otariidae: {otariidae_count}')
print(f'canidae: {canidae_count}')
print(f'mustelidae: {mustelidae_count}')
print(f'physeteridae: {physeteridae_count}')
print(f'balaenidae: {balaenidae_count}')
print(f'phocoenidae: {phocoenidae_count}')
print(f'rhinocerotidae: {rhinocerotidae_count}')
print(f'balaenopteridae: {balaenopteridae_count}')
print(f'monodontidae: {monodontidae_count}')
print(f'pontoporiidae: {pontoporiidae_count}')
print(f'lipotidae: {lipotidae_count}')
print(f'camelidae: {camelidae_count}')
print(f'ziphiidae: {ziphiidae_count}')
print(f'hippopotamidae: {hippopotamidae_count}')
print(f'vespertilionidae: {vespertilionidae_count}')
print(f'suidae: {suidae_count}')
print(f'pteropodidae: {pteropodidae_count}')
print(f'equidae: {equidae_count}')
print(f'cervidae: {cervidae_count}')
print(f'tapiridae: {tapiridae_count}')
print(f'bovidae: {bovidae_count}')
print(f'rhinolopidae: {rhinolopidae_count}')
print(f'moschidae: {moschidae_count}')
print(f'indriidae: {indriidae_count}')
print(f'leporidae: {leporidae_count}')
print(f'cheirogaleidae: {cheirogaleidae_count}')
print(f'cercopithecidae: {cercopithecidae_count}')
print(f'dasypodidae: {dasypodidae_count}')
print(f'hominidae: {hominidae_count}')
print(f'sciuridae: {sciuridae_count}')
print(f'callitrichidae: {callitrichidae_count}')


total_count_6 = (
    felidae_count + hyaenidae_count + delphinidae_count + ursidae_count + herpestidae_count +
    odobenidae_count + phocidae_count + otariidae_count + canidae_count + mustelidae_count +
    physeteridae_count + balaenidae_count + phocoenidae_count + rhinocerotidae_count +
    balaenopteridae_count + monodontidae_count + pontoporiidae_count + lipotidae_count +
    camelidae_count + ziphiidae_count + hippopotamidae_count + vespertilionidae_count +
    suidae_count + pteropodidae_count + equidae_count + cervidae_count + tapiridae_count +
    bovidae_count + rhinolopidae_count + moschidae_count + indriidae_count + leporidae_count +
    cheirogaleidae_count + cercopithecidae_count + dasypodidae_count + hominidae_count + sciuridae_count + callitrichidae_count
)

print(f'Total de contagens: {total_count_6}')

for seq_record in SeqIO.parse('seqdump (20).txt', 'fasta'):
    s_20 = ' '
    seqdump_20 = seq_record.id
    seqdump_20_seq = seq_record.seq
    description_20 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_20:
        description_20 = description_20.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_20 = description_20.split()[1:3]
    names_seq_20 = s_20.join(names_20)
    
    if names_seq_20 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_20)
        lista.append(f'>{names_seq_20}')
        lista.append(str(seqdump_20_seq))
for seq_record in SeqIO.parse('seqdump (21).txt', 'fasta'):
    s_21 = ' '
    seqdump_21 = seq_record.id
    seqdump_21_seq = seq_record.seq
    description_21 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_21:
        description_21 = description_21.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_21 = description_21.split()[1:3]
    names_seq_21 = s_21.join(names_21)
    
    if names_seq_21 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_21)
        lista.append(f'>{names_seq_21}')
        lista.append(str(seqdump_21_seq))
for seq_record in SeqIO.parse('seqdump (22).txt', 'fasta'):
    s_22 = ' '
    seqdump_22 = seq_record.id
    seqdump_22_seq = seq_record.seq
    description_22 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_22:
        description_22 = description_22.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_22 = description_22.split()[1:3]
    names_seq_22 = s_22.join(names_22)
    
    if names_seq_22 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_22)
        lista.append(f'>{names_seq_22}') # Macaca mulatta e Macaca fuscata MAIS UM
        lista.append(str(seqdump_22_seq))
    if names_seq_22 in felidae:
            felidae_count += 1
    if names_seq_22 in hyaenidae:
            hyaenidae_count += 1
    if names_seq_22 in delphinidae:
            delphinidae_count += 1
    if names_seq_22 in ursidae:
            ursidae_count += 1
    if names_seq_22 in herpestidae:
        herpestidae_count += 1
    if names_seq_22 in odobenidae:
        odobenidae_count += 1
    if names_seq_22 in phocidae:
        phocidae_count += 1
    if names_seq_22 in otariidae:
        otariidae_count += 1
    if names_seq_22 in canidae:
        canidae_count += 1
    if names_seq_22 in mustelidae:
        mustelidae_count += 1
    if names_seq_22 in physeteridae:
        physeteridae_count += 1
    if names_seq_22 in balaenidae:
        balaenidae_count += 1
    if names_seq_22 in phocoenidae:
        phocoenidae_count += 1
    if names_seq_22 in rhinocerotidae:
        rhinocerotidae_count += 1
    if names_seq_22 in balaenopteridae:
        balaenopteridae_count += 1
    if names_seq_22 in monodontidae:
        monodontidae_count += 1
    if names_seq_22 in pontoporiidae:
        pontoporiidae_count += 1
    if names_seq_22 in lipotidae:
        lipotidae_count += 1
    if names_seq_22 in camelidae:
        camelidae_count += 1
    if names_seq_22 in ziphiidae:
        ziphiidae_count += 1
    if names_seq_22 in hippopotamidae:
        hippopotamidae_count += 1
    if names_seq_22 in vespertilionidae:
        vespertilionidae_count += 1
    if names_seq_22 in suidae:
        suidae_count += 1
    if names_seq_22 in pteropodidae:
        pteropodidae_count += 1
    if names_seq_22 in equidae:
        equidae_count += 1
    if names_seq_22 in cervidae:
        cervidae_count += 1
    if names_seq_22 in tapiridae:
        tapiridae_count += 1
    if names_seq_22 in bovidae:
        bovidae_count += 1
    if names_seq_22 in rhinolopidae:
        rhinolopidae_count += 1
    if names_seq_22 in moschidae:
        moschidae_count += 1
    if names_seq_22 in indriidae:
        indriidae_count += 1
    if names_seq_22 in leporidae:
        leporidae_count += 1
    if names_seq_22 in cheirogaleidae:
        cheirogaleidae_count += 1
    if names_seq_22 in cercopithecidae:
        cercopithecidae_count += 1
    if names_seq_22 in dasypodidae:
        dasypodidae_count += 1
    if names_seq_22 in hominidae:
        hominidae_count += 1
    if names_seq_22 in sciuridae:
        sciuridae_count += 1
    if names_seq_22 in callitrichidae:
        callitrichidae_count += 1
print(f'felidae: {felidae_count}')
print(f'hyaenidae: {hyaenidae_count}')
print(f'delphinidae: {delphinidae_count}')
print(f'ursidae: {ursidae_count}')
print(f'herpestidae: {herpestidae_count}')
print(f'odobenidae:{odobenidae_count}')
print(f'phocidae: {phocidae_count}')
print(f'otariidae: {otariidae_count}')
print(f'canidae: {canidae_count}')
print(f'mustelidae: {mustelidae_count}')
print(f'physeteridae: {physeteridae_count}')
print(f'balaenidae: {balaenidae_count}')
print(f'phocoenidae: {phocoenidae_count}')
print(f'rhinocerotidae: {rhinocerotidae_count}')
print(f'balaenopteridae: {balaenopteridae_count}')
print(f'monodontidae: {monodontidae_count}')
print(f'pontoporiidae: {pontoporiidae_count}')
print(f'lipotidae: {lipotidae_count}')
print(f'camelidae: {camelidae_count}')
print(f'ziphiidae: {ziphiidae_count}')
print(f'hippopotamidae: {hippopotamidae_count}')
print(f'vespertilionidae: {vespertilionidae_count}')
print(f'suidae: {suidae_count}')
print(f'pteropodidae: {pteropodidae_count}')
print(f'equidae: {equidae_count}')
print(f'cervidae: {cervidae_count}')
print(f'tapiridae: {tapiridae_count}')
print(f'bovidae: {bovidae_count}')
print(f'rhinolopidae: {rhinolopidae_count}')
print(f'moschidae: {moschidae_count}')
print(f'indriidae: {indriidae_count}')
print(f'leporidae: {leporidae_count}')
print(f'cheirogaleidae: {cheirogaleidae_count}')
print(f'cercopithecidae: {cercopithecidae_count}')
print(f'dasypodidae: {dasypodidae_count}')
print(f'hominidae: {hominidae_count}')
print(f'sciuridae: {sciuridae_count}')
print(f'callitrichidae: {callitrichidae_count}')


total_count_5 = (
    felidae_count + hyaenidae_count + delphinidae_count + ursidae_count + herpestidae_count +
    odobenidae_count + phocidae_count + otariidae_count + canidae_count + mustelidae_count +
    physeteridae_count + balaenidae_count + phocoenidae_count + rhinocerotidae_count +
    balaenopteridae_count + monodontidae_count + pontoporiidae_count + lipotidae_count +
    camelidae_count + ziphiidae_count + hippopotamidae_count + vespertilionidae_count +
    suidae_count + pteropodidae_count + equidae_count + cervidae_count + tapiridae_count +
    bovidae_count + rhinolopidae_count + moschidae_count + indriidae_count + leporidae_count +
    cheirogaleidae_count + cercopithecidae_count + dasypodidae_count + hominidae_count + sciuridae_count + callitrichidae_count
)

print(f'Total de contagens: {total_count_5}')

for seq_record in SeqIO.parse('seqdump (23).txt', 'fasta'):
    s_23 = ' '
    seqdump_23 = seq_record.id
    seqdump_23_seq = seq_record.seq
    description_23 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_23:
        description_23 = description_23.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_23 = description_23.split()[1:3]
    names_seq_23 = s_23.join(names_23)
    
    if names_seq_23 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_23)
        lista.append(f'{names_seq_23}')
for seq_record in SeqIO.parse('seqdump (24).txt', 'fasta'):
    s_24 = ' '
    seqdump_24 = seq_record.id
    seqdump_24_seq = seq_record.seq
    description_24 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_24:
        description_24 = description_24.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_24 = description_24.split()[1:3]
    names_seq_24 = s_24.join(names_24)
    
    if names_seq_24 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_24)
        lista.append(f'{names_seq_24}')
for seq_record in SeqIO.parse('seqdump (25).txt', 'fasta'):
    s_25 = ' '
    seqdump_25 = seq_record.id
    seqdump_25_seq = seq_record.seq
    description_25 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_25:
        description_25 = description_25.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_25 = description_25.split()[1:3]
    names_seq_25 = s_25.join(names_25)
    
    if names_seq_25 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_25)
        lista.append(f'{names_seq_25}')
for seq_record in SeqIO.parse('seqdump (26).txt', 'fasta'):
    s_26 = ' '
    seqdump_26 = seq_record.id
    seqdump_26_seq = seq_record.seq
    description_26 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_26:
        description_26 = description_26.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_26 = description_26.split()[1:3]
    names_seq_26 = s_26.join(names_26)
    
    if names_seq_26 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_26)
        lista.append(f'{names_seq_26}')
for seq_record in SeqIO.parse('seqdump (27).txt', 'fasta'):
    s_27 = ' '
    seqdump_27 = seq_record.id
    seqdump_27_seq = seq_record.seq
    description_27 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_27:
        description_27 = description_27.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_27 = description_27.split()[1:3]
    names_seq_27 = s_27.join(names_27)
    
    if names_seq_27 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_27)
        lista.append(f'{names_seq_27}')
for seq_record in SeqIO.parse('seqdump (28).txt', 'fasta'):
    s_28 = ' '
    seqdump_28 = seq_record.id
    seqdump_28_seq = seq_record.seq
    description_28 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_28:
        description_28 = description_28.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_28 = description_28.split()[1:3]
    names_seq_28 = s_28.join(names_28)
    
    if names_seq_28 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_28)
        lista.append(f'{names_seq_28}')
for seq_record in SeqIO.parse('seqdump (29).txt', 'fasta'):
    s_29 = ' '
    seqdump_29 = seq_record.id
    seqdump_29_seq = seq_record.seq
    description_29 = seq_record.description
    
    # Removendo "PREDICTED:" se estiver presente
    if 'PREDICTED:' in description_29:
        description_29 = description_29.replace('PREDICTED:', '')

    # Capturando o nome da espécie a partir da descrição modificada
    names_29 = description_29.split()[1:3]
    names_seq_29 = s_29.join(names_29)
    
    if names_seq_29 in unique_ids:
        continue  # Ignorar duplicados
    else:
        unique_ids.add(names_seq_29)
        lista.append(f'{names_seq_29}')
with open('nomes_mc1r_seq.fasta', 'w') as arquivo:
    for item in lista:
        arquivo.write(item + '\n')

warnings.simplefilter('ignore', DeprecationWarning)
aligner = Align.PairwiseAligner()

fasta_file = "/home/caiog/PIBIC/mc1r_felinae_blast/nomes_mc1r_seq.fasta"

# Comando para rodar o MAFFT
mafft_cline = MafftCommandline(input=fasta_file)

# Executando o comando e capturando a saída
stdout, stderr = mafft_cline()

with open("nomes_mafft.fasta", "w") as aligned_file:
    aligned_file.write(stdout)

# Carregar o arquivo de alinhamento
alignment = AlignIO.read("nomes_mafft.fasta", "fasta")

# Aparar o alinhamento (mantendo as colunas de 5852 a 6614)
trimmed_alignment = alignment[:, 5852:6614]

# Definir o padrão de interesse (exatamente 24 '-')
pattern = '-' * 24

# Verificar cada sequência no alinhamento aparado
for seq_index, record in enumerate(trimmed_alignment):
    seq_str = str(record.seq)
    if pattern in seq_str and seq_str.count('-') == 24:
        print(f"Padrão de 24 '-' encontrado na sequência {seq_index + 1} {record.description}")

