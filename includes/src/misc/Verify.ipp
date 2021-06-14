/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTION: VERIFY_OPEN
          This function verifies the status of a inputfile strem.
          If the file is not open for reading, it gives ERROR and abort.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void mpm::misc::VERIFY_OPEN(std::ifstream& file, std::string& fileName) {
    if (!file.is_open()) {
        std::cerr << "ERROR: in opening file " << fileName << "\n";
        abort();
    }
}
