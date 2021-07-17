// InputParser
//
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
//
// ABOUT
// - Parses JSON-like input files with parameters and organizes the input
//   into a ParameterPack (see class comments, below)


#ifndef INPUT_PARSER_H
#define INPUT_PARSER_H

#include <fstream>
#include <exception>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "OrderParameters/tools/StringTools.h"


// ParameterPack: 
// - data structure which organizes parsed input
//   - For simplicity, I recommend using the provided get functions to access input
//     rather than accessing the data directly
//   - However, all member variables are public for convenience
// - Stores all values (whether or not they're numbers) as strings, for generality
// - Data members
//   - values: key-value pairs
//   - vectors: keys which map to vectors of strings
//   - parameter_packs: sub-ParameterPacks which provide more coarse-grained groupings
//     of parameters
//     - Useful for passing only *parts* of an input file to the constructor that needs them
// - Member functions
//   - findValues, findVectors, findParameterPacks  (plural)
//     - Return a vector with ptrs to all of the matches from the corresponding multimap
//   - findValue, findVector, findParameterPack  (singular)
//     - Expect a unique entry with the given key, and throw if there are duplicates
//     - Throw if the KeyType is Required, and no matches are found
//
// TODO:
// - Functions that return by value for cleaner semantics?
class ParameterPack
{
 public:
  ParameterPack() = default;

  // Creates an empty pack with the given name
  ParameterPack(const std::string& pack_name);
 
  // Indicates whether the given key must be present
  enum class KeyType {
    Required,
    Optional
  };


  //----- Adding Entries ----//

  // Inserts a copy of 'value' under entry 'key' and returns a reference to the copy
  std::string& insert(
    const std::string& key,
    const std::string& value
  );
  
  // Inserts a copy of 'vec' under entry 'key' and returns a reference to the copy
  std::vector<std::string>& insert(
    const std::string& key,
    const std::vector<std::string>& vec
  );

  // Inserts a copy of 'pack' under entry 'key' and returns a reference to the copy
  ParameterPack& insert(
    const std::string& key,
    const ParameterPack& pack
  );


  //----- Removing Entries -----//

  // Erases all values registered under 'key'
  void eraseValues(const std::string& key) {
    values_.erase(key);
  }


  //----- Input Conversion Functions -----//

  // Return 'true' if a match was found and conversion was successful

  // Reads a unique number from the ParameterPack and performs the
  // necessary conversions
  template<typename T>
  bool readNumber(const std::string& key, const KeyType& key_type, T& value) const {
    const std::string* value_ptr = findValue(key, key_type);
    if ( value_ptr != nullptr ) {
      value = StringTools::stringToValue<T>(*value_ptr);
      return true;
    }
    else {
      return false;
    }
  }

  // Reads a unique flag (a bool) from the ParameterPack and performs the
  // necessary conversions
  bool readFlag(const std::string& key, const KeyType& key_type, bool& flag) const {
    const std::string* value_ptr = findValue(key, key_type);
    if ( value_ptr != nullptr ) {
      flag = StringTools::stringToBool(*value_ptr);
      return true;
    }
    else {
      return false;
    }
  }

  bool readBool(const std::string& key, const KeyType& key_type, bool& flag) {
    return readFlag(key, key_type, flag);
  }

  // Reads a unique string from the ParameterPack
  bool readString(const std::string& key, const KeyType& key_type, std::string& str) const {
    const std::string* value_ptr = findValue(key, key_type);
    if ( value_ptr != nullptr ) {
      // No conversion necessary
      str = *value_ptr;
      return true;
    }
    else {
      return false;
    }
  }

  // Reads a unique array from the ParameterPack, converting from std::string
  // to the appropriate type
  template<typename T, std::size_t Dim>
  bool readArray(const std::string& key, const KeyType& key_type, std::array<T,Dim>& arr) const {
    const auto* vec_ptr = findVector(key, key_type);
    if ( vec_ptr != nullptr ) {
      arr = StringTools::stringsToArray<T,Dim>(*vec_ptr);
      return true;
    }
    else {
      return false;
    }
  }

  // Reads a unique vector from the ParameterPack, converting from std::string
  // to the appropriate type
  template<typename T>
  bool readVector(const std::string& key, const KeyType& key_type, std::vector<T>& vec) const {
    const auto* vec_ptr = findVector(key, key_type);
    if ( vec_ptr != nullptr ) {
      vec = StringTools::stringsToVector<T>(*vec_ptr);
      return true;
    }
    else {
      return false;
    }
  }


  //----- Searching Functions -----//

  // Find unique entries
  // - These functions throw if more than one entry is found for the given key
  // - Note that these functions call their non-unique-entry versions, which throw 
  //   if KeyType is "Required" and the key is not found
  const std::string* findValue(
    const std::string& key,
    const KeyType&     key_type
  ) const;

  const std::vector<std::string>* findVector(
    const std::string& key,
    const KeyType&     key_type
  ) const;

  const ParameterPack* findParameterPack(
    const std::string& key,
    const KeyType&     key_type
  ) const; 

  // Return by const ref (instead of const ptr) is only safe if the entries
  // are required
  const std::string& findRequiredValue(const std::string& key) const {
    return *findValue(key, KeyType::Required);
  }

  const std::vector<std::string>& findRequiredVector(const std::string& key) const {
    return *findVector(key, KeyType::Required);
  }

  const ParameterPack& findRequiredParameterPack(const std::string& key) const {
    return *findParameterPack(key, KeyType::Required);
  }

  // Find *all* entries matching the given key
  // - These functions return vectors of pointers to the matching entries
  //   - These pointers are guaranteed to be non-null
  std::vector<const std::string*> findValues(
    const std::string& key,
    const KeyType&     key_type
  ) const;
  std::vector<const std::vector<std::string>*> findVectors(
    const std::string& key,
    const KeyType&     key_type
  ) const;
  std::vector<const ParameterPack*> findParameterPacks(
    const std::string& key,
    const KeyType&     key_type
  ) const;


  //----- Misc. -----//

  void clear() {
    this->values_.clear();
    this->vectors_.clear();
    this->parameter_packs_.clear();
  };

  // Print all ParameterPack contents
  std::ostream& print(
    std::ostream& os,
    const std::string& indent = ""
  ) const;

  const std::string& getName() const {
    return name_;
  }


  //----- Exceptions -----//

  // A mandatory key was not found
  class KeyNotFoundException : public std::exception {
   public:
    KeyNotFoundException(const std::string& key, const std::string& context) {
      std::stringstream err_ss;
      err_ss << "Mandatory key \"" << key << "\" was not found in context \"" << context << "\"\n";
      message_ = err_ss.str();
    };

    const char* what() const noexcept override { return message_.c_str(); };

   private:
    std::string message_;
  };

  // A key is not unique
  class DuplicateKeyException : public std::exception {
   public:
    DuplicateKeyException(const std::string& key, const int num_duplicates,
                          const std::string& context) {
      std::stringstream err_ss;
      err_ss << "Unique key \"" << key << "\" appears " << num_duplicates 
             << " times in context \"" << context << "\"\n"
             << "  It is unclear which one should be used.\n";
      message_ = err_ss.str();
    };

    const char* what() const noexcept override { return message_.c_str(); };

   private:
    std::string message_;
  };

 private:
  // Name of this organization unit
  std::string name_;

  // Parameters
  std::multimap<std::string, std::string>              values_;
  std::multimap<std::string, std::vector<std::string>> vectors_;
  std::multimap<std::string, ParameterPack>            parameter_packs_;
};


// Reads whitespace-delimited tokens from a std::iostream, ignoring comments
// TODO:
// - Custom comment characters
// - Allow omission of whitespace between certain characters to make input
//   even more forgiving
class TokenStream {
 public:
  // Create a TokenStream from a std::iostream
  // - Note that TokenStream does not claim ownership of the iostream,
  //   but changes the iostream's state as the TokenStream reads tokens from it
  TokenStream(std::iostream& io_stream) 
   : io_stream_(io_stream) 
  {}

  enum class Status : int {
    Success,
    EndOfStream,    // e.g. EOF
    ClosingBrace,   // right curly brace, }
    ClosingBracket, // right square bracket, ]
    Failure
  };

  // Get the next token from the TokenStream, discarding any comments
  // along the way
  TokenStream::Status getNextToken(
    std::string& token
  );

  std::string getCurrentLine() const {
    return line_stream_.str();
  }

 private:
  std::iostream&    io_stream_;
  std::stringstream line_stream_;

  std::string comment_chars_ = "#";
};



// Reads and parses an OrderParameters input file
class InputParser
{
 public:
  InputParser() = default;

  void parseFile(
    const std::string file,
    ParameterPack& parameter_pack
  ) const;

  ParameterPack parseFile(const std::string file) const {
    ParameterPack pack;
    parseFile(file, pack);
    return pack;
  }

 private:
  // Characters which indicate the beginning of comments
  std::string comment_chars_ = "#";


  // Read the next ParameterPack entry from the TokenStream
  // and add it to the given ParameterPack
  TokenStream::Status parseNextEntry(
    TokenStream&   token_stream,
    ParameterPack& parameter_pack
  ) const;

  // Parses a vector of strings from the TokenStream
  // - Assumes that the opening bracket, [,  has already been parsed
  //   -  TODO place the bracket back into the line_stream before calling?
  // - Reads until encountering the next closing bracket, ]
  TokenStream::Status parseVector(
    TokenStream& token_stream,
    std::vector<std::string>& vec
  ) const;

  // Parses a ParameterPack from the TokenStream
  // - Assumes that the opening brace, {,  has already been parsed
  //   -  TODO place the brace back into the line_stream before calling?
  // - Reads until encountering the closing brace, }, in its scope
  TokenStream::Status parseParameterPack(
    TokenStream&   token_stream,
    ParameterPack& parameter_pack
  ) const;


  //----- Stream Manipulation -----//

  // Puts the remaining characters in the stringstream back into the iostream
  // used to generate it
  // - "line_stream" is a partially-read stringstream that contains a line of text in 
  //   which was generated from "stream" using getline()
  // - Assumes:
  //   - that getline() has not been called again since line_stream was created
  //   - that the newline character ('\n') is the delimiter
  void putBackRestOfLine(std::iostream& io_stream, std::stringstream& line_stream) const {
    io_stream.putback('\n');
    for ( int i = line_stream.str().length() - 1; i >= line_stream.tellg(); --i ) {
      io_stream.putback( line_stream.str()[i] );
    }
  }

  // Unget the indicated number of characters obtained from the stream
  void ungetCharacters(std::iostream& stream, const int num_characters) const {
    for ( int i = num_characters; i >= 0; --i ) { 
      stream.unget(); 
    }
  };

  // Same as putBackRestOfLine(), but uses unget
  void ungetRestOfLine(std::iostream& io_stream, std::stringstream& line_stream) const {
    ungetCharacters( io_stream, line_stream.str().length() - line_stream.tellg() );
  }


  //----- Exceptions -----//

  // Incomplete entry in a ParameterPack
  class IncompleteEntryException : public std::exception {
   public:
    IncompleteEntryException(const std::string& key, const TokenStream& token_stream) {
      std::stringstream err_ss;
      err_ss << "key \"" << key << "\" has incomplete declaration\n" 
             << "  current line: " << token_stream.getCurrentLine() << '\n'
             << DebugTools::stacktrace() << '\n';
      message_ = err_ss.str();
    };

    const char* what() const noexcept override { return message_.c_str(); };

   private:
    std::string message_;
  };

  // Missing delimiter in a ParameterPack entry
  class MissingDelimiterException : public std::exception {
   public:
    MissingDelimiterException(const std::string& key, const TokenStream& token_stream) {
      std::stringstream err_ss;
      err_ss << "key \"" << key << "\" is missing its \"=\" delimiter\n"
             << "  current line: " << token_stream.getCurrentLine() << '\n'
             << DebugTools::stacktrace() << '\n';
      message_ = err_ss.str();
    };

    const char* what() const noexcept override { return message_.c_str(); };

   private:
    std::string message_;
  };

  // Vector missing its closing bracket
  class MissingBracketException : public std::exception {
   public:
    MissingBracketException(const std::string& name, const TokenStream& token_stream) {
      std::stringstream err_ss;
      err_ss << "Vector \"" << name << "\" is missing its closing bracket\n"
             << "  current line: " << token_stream.getCurrentLine() << '\n'
             << DebugTools::stacktrace() << '\n';
      message_ = err_ss.str();
    };

    const char* what() const noexcept override { return message_.c_str(); };

   private:
    std::string message_;
  };

  // ParameterPack missing its closing brace
  class MissingBraceException : public std::exception {
   public:
    MissingBraceException(const std::string& name, const TokenStream& token_stream) {
      std::stringstream err_ss;
      err_ss << "ParamterPack \"" << name << "\" is missing its closing brace\n"
             << "  current line: " << token_stream.getCurrentLine() << '\n'
             << DebugTools::stacktrace() << '\n';
      message_ = err_ss.str();
    };

    const char* what() const noexcept override { return message_.c_str(); };

   private:
    std::string message_;
  };
};

#endif // ifndef INPUT_PARSER_H