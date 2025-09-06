// CMSC 341 - Spring 25 - Project 4
#include "dnadb.h"
DnaDb::DnaDb(int size, hash_fn hash, prob_t probing = DEFPOLCY){
    if (isPrime(size) && size >= MINPRIME && size <= MAXPRIME){
        m_currentCap = size;
    }else{
        m_currentCap = findNextPrime(size);
    }

    m_currentTable = new DNA*[m_currentCap];
    for(int i = 0; i < m_currentCap; i++){
        m_currentTable[i] = nullptr;
    }
    m_currentSize = 0;
    m_hash = hash;
    m_currProbing = probing;
    m_newPolicy = probing;
    m_currNumDeleted = 0;

    //Setting old variables to the new table since a transfer would have the same cap
    m_oldCap = m_currentCap;
    m_oldTable = nullptr;
    m_oldSize = 0;
    m_oldProbing = probing;
    m_oldNumDeleted = 0;

    m_transferIndex = 0;
}

DnaDb::~DnaDb(){

    if(m_oldTable != nullptr){
        for(int i = 0; i < m_oldCap; i++){
            if(m_oldTable[i] != nullptr){
                delete m_oldTable[i];
                m_oldTable[i] = nullptr;
            }
        }
        delete[] m_oldTable;
        m_oldTable = nullptr;
    }

    if(m_currentTable != nullptr){
        for(int i = 0; i < m_currentCap; i++){
            if(m_currentTable[i] != nullptr){
                delete m_currentTable[i];
                m_currentTable[i] = nullptr;
            }
        }
        delete[] m_currentTable;
        m_currentTable = nullptr;
    }
}

void DnaDb::changeProbPolicy(prob_t policy){
    m_newPolicy = policy;
}

bool DnaDb::insert(DNA dna){
    int insertIndex;
    int i = 0;
    bool inserted = false;
    insertIndex = (m_hash(dna.getSequence()) % m_currentCap);
    if(getDNA(dna.getSequence(), dna.getLocId()).getUsed() == true || dna.getLocId() < MINLOCID || dna.getLocId() > MAXLOCID){
        return false;
    }
    //Loops till the node is inserted since we already duplicate checked above
    while(!inserted){
        i++;
        if(m_currentTable[insertIndex] == nullptr){
            m_currentTable[insertIndex] = new DNA(dna.getSequence(), dna.getLocId(), true);
            inserted = true;
            m_currentSize++;
        }else if(m_currentTable[insertIndex]->getUsed() == false){
            m_currentTable[insertIndex]->setLocID(dna.getLocId());
            m_currentTable[insertIndex]->setSequence(dna.getSequence());
            m_currentTable[insertIndex]->setUsed(true);
            inserted = true;
            m_currNumDeleted--;
        }else{
            switch (m_currProbing)
            {
            case QUADRATIC:
                insertIndex = ((m_hash(dna.getSequence()) + (i * 2)) % m_currentCap); 
                break;
            case LINEAR:
                insertIndex = ((m_hash(dna.getSequence()) + i) % m_currentCap);
                break;
            case DOUBLEHASH:
                insertIndex = ((m_hash(dna.getSequence()) + (i * (11 - ((m_hash(dna.getSequence()) % 11))))) % m_currentCap);
                break;
            default:
                break;
            }
        }
    }
    //rehash checker
    if(!inserted){return false;}

    if(lambda() > MAXLOAD || m_oldTable != nullptr){
        rehash();
    }

    return inserted;


}

void DnaDb::rehash(){
    //If first rehash move curr table to old table
    if(m_oldTable == nullptr){
        m_oldTable = m_currentTable;
        m_oldCap = m_currentCap;
        m_oldSize = m_currentSize;
        m_oldProbing = m_currProbing;
        m_oldNumDeleted = m_currNumDeleted;

        m_currentCap = findNextPrime((m_currentSize - m_currNumDeleted) * 4);
        m_currentTable = new DNA*[m_currentCap];
        m_currentSize = 0;
        m_currNumDeleted = 0;
        m_currProbing = m_newPolicy;
    
        for(int i = 0; i < m_currentCap; i++){
            m_currentTable[i] = nullptr;
        }

    }
    bool rehashDone = false;
    int quarterHash = floor(m_transferIndex + (m_oldCap / 4));
    if(quarterHash > (3 * m_oldCap) / 4){
        quarterHash = m_oldCap;
        rehashDone = true;
    }
    //trasnfer over 1/4th the data at a time
    for(; m_transferIndex < quarterHash; m_transferIndex++){
        if(m_oldTable[m_transferIndex] != nullptr && m_oldTable[m_transferIndex]->getUsed()){
            m_oldTable[m_transferIndex]->setUsed(false);
            insert(*m_oldTable[m_transferIndex]);
        }
    }
    //once rehashing is complete clear the old table
    if(rehashDone){
        for(int i = 0; i < m_oldCap; i++){
            delete m_oldTable[i];
            m_oldTable[i] = nullptr;
        }
        delete[] m_oldTable;
        m_oldTable = nullptr;
    }
}

bool DnaDb::remove(DNA dna){
    int insertIndexCurr, insertIndexOld;
    int i = 0;
    bool dnaFound = false;
    insertIndexCurr = (m_hash(dna.getSequence()) % m_currentCap);
    
    // Search in current table till nullptr or the dna
    while(m_currentTable[insertIndexCurr] != nullptr && !dnaFound){
        if(m_currentTable[insertIndexCurr]->getSequence() == dna.getSequence() && 
            m_currentTable[insertIndexCurr]->getLocId() == dna.getLocId()){
            m_currentTable[insertIndexCurr]->setUsed(false);
            m_currNumDeleted ++;
            dnaFound = true;
            break;
        }
        i++;
        switch (m_currProbing)
        {
        case QUADRATIC:
            insertIndexCurr = ((m_hash(dna.getSequence()) + (i * 2)) % m_currentCap); 
            break;
        case LINEAR:
            insertIndexCurr = ((m_hash(dna.getSequence()) + i) % m_currentCap);
            break;
        case DOUBLEHASH:
            insertIndexCurr = ((m_hash(dna.getSequence()) + (i * (11 - ((m_hash(dna.getSequence()) % 11))))) % m_currentCap);
            break;
        default:
            break;
        }
    }
    //check if old table exist then search
    if(m_oldTable != nullptr && m_oldCap > 0){
        insertIndexOld = (m_hash(dna.getSequence()) % m_oldCap);
        i = 0;
        //loop till either node is found or an empty slot is found
        while(m_oldTable[insertIndexOld] != nullptr){
            if(m_oldTable[insertIndexOld]->getSequence() == dna.getSequence() && 
            m_oldTable[insertIndexOld]->getLocId() == dna.getLocId()){
                m_oldTable[insertIndexOld]->setUsed(false);
                dnaFound = true;
                m_oldNumDeleted ++;
                break;
            }
            i++;
            //probing
            switch (m_oldProbing)
            {
            case QUADRATIC:
                insertIndexOld = ((m_hash(dna.getSequence()) + (i * i)) % m_oldCap); 
                break;
            case LINEAR:
                insertIndexOld = ((m_hash(dna.getSequence()) + i) % m_oldCap);
                break;
            case DOUBLEHASH:
                insertIndexOld = ((m_hash(dna.getSequence()) + (i * (11 - ((m_hash(dna.getSequence()) % 11))))) % m_oldCap);
                break;
            default:
                break;
            }
        }
    }
    //if not inserted return false else check to see if rehasing must me done
    if(!dnaFound){
        return false;
    }else if(deletedRatio() > MAXDELETED || m_oldTable != nullptr){
        rehash();
    }
    return dnaFound;
}

const DNA DnaDb::getDNA(string sequence, int location) const{
    int insertIndexCurr;
    int i = 0;
    insertIndexCurr = (m_hash(sequence) % m_currentCap);
    DNA empty = DNA();
    if(location < MINLOCID || location > MAXLOCID){
        return empty;
    }
    // Search in current table
    while(m_currentTable[insertIndexCurr] != nullptr){
        if(m_currentTable[insertIndexCurr]->getSequence() == sequence && 
           m_currentTable[insertIndexCurr]->getLocId() == location){
            return *m_currentTable[insertIndexCurr];
        }
        i++;
        switch (m_currProbing)
        {
        case QUADRATIC:
            insertIndexCurr = ((m_hash(sequence) + (i * i)) % m_currentCap); 
            break;
        case LINEAR:
            insertIndexCurr = ((m_hash(sequence) + i) % m_currentCap);
            break;
        case DOUBLEHASH:
            insertIndexCurr = ((m_hash(sequence) + (i * (11 - ((m_hash(sequence) % 11))))) % m_currentCap);
            break;
        default:
            break;
        }
    }
    
    // Search in old table only if it exists
    if(m_oldTable != nullptr && m_oldCap > 0){
        int insertIndexOld = (m_hash(sequence) % m_oldCap);
        i = 0;
        while(m_oldTable[insertIndexOld] != nullptr){
            if(m_oldTable[insertIndexOld]->getSequence() == sequence && 
               m_oldTable[insertIndexOld]->getLocId() == location){
                return *m_oldTable[insertIndexOld];
            }
            i++;
            switch (m_oldProbing)
            {
            case QUADRATIC:
                insertIndexOld = ((m_hash(sequence) + (i * i)) % m_oldCap); 
                break;
            case LINEAR:
                insertIndexOld = ((m_hash(sequence) + i) % m_oldCap);
                break;
            case DOUBLEHASH:
                insertIndexOld = ((m_hash(sequence) + (i * (11 - ((m_hash(sequence) % 11))))) % m_oldCap);
                break;
            default:
                break;
            }
        }
    }
    return empty; 
}

bool DnaDb::updateLocId(DNA dna, int location){
    int insertIndexCurr, insertIndexOld;
    int i = 0;
    insertIndexCurr = (m_hash(dna.getSequence()) % m_currentCap);

    
    // Search in current table
    while(m_currentTable[insertIndexCurr] != nullptr){
        if(m_currentTable[insertIndexCurr]->getSequence() == dna.getSequence() && 
           m_currentTable[insertIndexCurr]->getLocId() == dna.getLocId()){
            m_currentTable[insertIndexCurr]->setLocID(location);
            return true;
        }
        i++;
        switch (m_currProbing)
        {
        case QUADRATIC:
            insertIndexCurr = ((m_hash(dna.getSequence()) + (i * 2)) % m_currentCap); 
            break;
        case LINEAR:
            insertIndexCurr = ((m_hash(dna.getSequence()) + i) % m_currentCap);
            break;
        case DOUBLEHASH:
            insertIndexCurr = ((m_hash(dna.getSequence()) + (i * (11 - ((m_hash(dna.getSequence()) % 11))))) % m_currentCap);
            break;
        default:
            break;
        }
    }
    
    // Search in old table
    if(m_oldTable != nullptr && m_oldCap > 0){
        insertIndexOld = (m_hash(dna.getSequence()) % m_oldCap);
        i = 0;
        while(m_oldTable[insertIndexOld] != nullptr){
            if(m_oldTable[insertIndexOld]->getSequence() == dna.getSequence() && 
            m_oldTable[insertIndexOld]->getLocId() == dna.getLocId()){
                m_oldTable[insertIndexOld]->setLocID(location);
                return true;
            }
            i++;
            switch (m_oldProbing)
            {
            case QUADRATIC:
                insertIndexOld = ((m_hash(dna.getSequence()) + (i * i)) % m_oldCap); 
                break;
            case LINEAR:
                insertIndexOld = ((m_hash(dna.getSequence()) + i) % m_oldCap);
                break;
            case DOUBLEHASH:
                insertIndexOld = ((m_hash(dna.getSequence()) + (i * (11 - ((m_hash(dna.getSequence()) % 11))))) % m_oldCap);
                break;
            default:
                break;
            }
        }
    }
    return false; 
}

float DnaDb::lambda() const {
      return (float(m_currentSize)/float(m_currentCap));
}

float DnaDb::deletedRatio() const {
    return (float(m_currNumDeleted)/float(m_currentSize));
}

void DnaDb::dump() const {
    cout << "Dump for the current table: " << endl;
    if (m_currentTable != nullptr)
        for (int i = 0; i < m_currentCap; i++) {
            cout << "[" << i << "] : " << m_currentTable[i] << endl;
        }
    cout << "Dump for the old table: " << endl;
    if (m_oldTable != nullptr)
        for (int i = 0; i < m_oldCap; i++) {
            cout << "[" << i << "] : " << m_oldTable[i] << endl;
        }
}

bool DnaDb::isPrime(int number){
    bool result = true;
    for (int i = 2; i <= number / 2; ++i) {
        if (number % i == 0) {
            result = false;
            break;
        }
    }
    return result;
}

int DnaDb::findNextPrime(int current){
    //we always stay within the range [MINPRIME-MAXPRIME]
    //the smallest prime starts at MINPRIME
    if (current < MINPRIME) current = MINPRIME-1;
    for (int i=current; i<MAXPRIME; i++) { 
        for (int j=2; j*j<=i; j++) {
            if (i % j == 0) 
                break;
            else if (j+1 > sqrt(i) && i != current) {
                return i;
            }
        }
    }
    //if a user tries to go over MAXPRIME
    return MAXPRIME;
}