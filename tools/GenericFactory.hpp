// GenericFactory
// AUTHOR: Sean M. Marks (https://github.com/seanmarks)
// - Implements a factory for registering and creating objects with a
//   shared base class, based on a key (typically a string)

#pragma once
#ifndef GENERIC_FACTORY_H
#define GENERIC_FACTORY_H

#include <exception>
#include <functional>
#include <iostream>
#include <memory>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>

#include "Assert.hpp"

// Sources:
//   Jim Hyslop and Herb Sutter 
//     http://www.drdobbs.com/conversations-abstract-factory-template/184403786
//
//   StackOverflow
//     https://stackoverflow.com/questions/24119808/derived-class-selection-based-on-string-c

template <class    BaseType,
          class    DerivedTypeKey,
          class... DerivedTypeArgs>
class GenericFactory
{
 public:
	// Typedefs
	using DerivedTypeConstructor = std::function< BaseType*(DerivedTypeArgs&&...) >;
	using Registry               = std::map< DerivedTypeKey, DerivedTypeConstructor >;

	// Returns a handle to the static GenericFactory for this type
	// - Uses the "construct on first use" idiom
	static GenericFactory& factory() { 
		static GenericFactory factory_instance;
		return factory_instance;
	};

	// Constructor and key operators (TODO?)
	//GenericFactory() {};
	//GenericFactory(const GenericFactory&);
	//GenericFactory& operator=(const GenericFactory&);

	// Registers methods to construct derived types
	// - Throws an exception if the mapping already exists (likely a programmer mistake)
	void registerCreateFunction(const DerivedTypeKey& key, DerivedTypeConstructor constructor_ptr) {
		// See whether the mapping already exists
		typename Registry::const_iterator registry_entry = registry_.find(key);
		FANCY_ASSERT( registry_entry == registry_.end(),
		              "GenericFactory key \"" << key << "\" already exists" );
		registry_.insert( std::make_pair(key, constructor_ptr) );
	}

	// Creates an instance of DerivedType by passing "input" to its constructor,
	// then returns a ptr (to BaseType!) to that instance
	BaseType* create(const DerivedTypeKey& key, DerivedTypeArgs&&... input) const {
    typename Registry::const_iterator registry_entry = registry_.find(key);
		FANCY_ASSERT( registry_entry != registry_.end(),
		              "GenericFactory key \"" << key << "\" not found" );
    return registry_entry->second( std::forward<DerivedTypeArgs>(input)... );
	}

	const Registry& get_registry() const {
		return registry_;
	}

 private:
	Registry registry_;
};


// Template-class method for registering <Base,Derived,Key,Input> tuples
template <class    BaseType, 
          class    DerivedType, 
          class    DerivedTypeKey,   // Key which maps to desired DerivedType
          class... DerivedTypeArgs>  // Types of inputs to DerivedType
class RegisterInFactory
{
 public:
	static_assert(std::is_base_of<BaseType, DerivedType>::value, "invalid usage" );

	// Constructor registers the map from "key" to a lambda that produces the derived-class object
	RegisterInFactory(const DerivedTypeKey& key) {
		GenericFactory<BaseType,DerivedTypeKey,DerivedTypeArgs...>::factory().registerCreateFunction(
			// Keyword
			key,
			// Lambda to produce the derived class instance
			[](DerivedTypeArgs&&... input){ return ( new DerivedType( std::forward<DerivedTypeArgs>(input)... ) ); }
		);
	}
};

#endif /* GENERIC_FACTORY_H */