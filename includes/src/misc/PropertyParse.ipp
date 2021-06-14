/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTION: PARSE_PROPERTIES
          This fuction first tokenize the given line using white space 
          as the delimiter.
          Then the value correspond to the each parameter is saved by its
          appropriate type.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void mpm::misc::PARSE_INPUT_PARAMETERS(std::string& iLine) {

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer; 
    boost::char_separator<char> separator(" ");
    tokenizer token(iLine, separator);
    boost::tokenizer<boost::char_separator<char> >::iterator parameter = token.begin();
    std::string par = *parameter;
    std::advance(parameter,1);

    if (par == "gravityFlag") {
        try {
            gravity = boost::lexical_cast<int>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse properties" << "\n";
            abort();
        }
    }
    if (par == "freeSurface") {
        try {
            freeSurface = boost::lexical_cast<int>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter freeSurface" << "\n";
            abort();
        }
    }

    if (par == "contact") {
        try {
            contact = boost::lexical_cast<int>(*parameter);
            }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter contact" << "\n";
            abort();
        }
    }

    if (par == "dt") {
        try {
            dt_ = boost::lexical_cast<double>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter dt" << "\n";
            abort();
        }
    }
    if (par == "totalTime") {
        try {
            total_time_ = boost::lexical_cast<double>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter totalTime" << "\n";
            abort();
        }
    }
    if (par == "loadStep") {
        try {
            load_step_ = boost::lexical_cast<double>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter loadStep" << "\n";
            abort();
        }
    }
    if (par == "timeFactor") {
        try {
            time_factor_ = boost::lexical_cast<unsigned>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter timeFactor" << "\n";
            abort();
        }
    }
    if (par == "scalarBeta") {
        try {
            scalar_beta_ = boost::lexical_cast<unsigned>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter scalarBeta" << "\n";
            abort();
        }
    }
    if (par == "dampingCoefficient") {
        try {
            damping_coefficient_ = boost::lexical_cast<double>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter scalarBeta" << "\n";
            abort();
        }
    }
    if (par == "constantPermeability") {
        try {
            permeability_ = boost::lexical_cast<double>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter scalarBeta" << "\n";
            abort();
        }
    }
    if (par == "NewmarkGamma") {
        try {
            Gamma = boost::lexical_cast<double>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter scalarBeta" << "\n";
            abort();
        }
    }
    if (par == "NewmarkBeta") {
        try {
            Beta = boost::lexical_cast<double>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter scalarBeta" << "\n";
            abort();
        }
    }
    if (par == "compressibility") {
        try {
            compressibility_ = boost::lexical_cast<double>(*parameter);
        }
        catch (const boost::bad_lexical_cast &) {
            std::cerr << "ERROR: failed to parse input parameter scalarBeta" << "\n";
            abort();
        }
    }
}


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTION: READ_PROPERTIES
          This fuction first tokenize the given line using white space 
          as the delimiter.
          Then the value and the property name is stored in a std::map.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void mpm::misc::READ_PROPERTIES(std::string& line) {

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer; 
    boost::char_separator<char> separator(" ");
    tokenizer token(line, separator);
    boost::tokenizer<boost::char_separator<char> >::iterator parameter = token.begin();
    std::string propName = *parameter;
    std::advance(parameter,1);

    double propValue = boost::lexical_cast<double>(*parameter);
    propertyList[propName] = propValue;

}
