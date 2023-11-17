#ifndef L1Trigger_L1THGCal_Binary_SA_h
#define L1Trigger_L1THGCal_Binary_SA_h

#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdint>

struct binary
{
  binary( const uint64_t& value ): value_( value ){}
  inline bool operator!= ( const binary& other ) const { return other.value_ != value_; }
  inline bool operator!= ( const uint64_t& other ) const { return other != value_; }
  inline bool operator== ( const binary& other ) const { return other.value_ == value_; }
  inline bool operator== ( const uint64_t& other ) const { return other == value_; }
  inline bool operator<= ( const binary& other ) const { return other.value_ <= value_; }
  inline bool operator<= ( const uint64_t& other ) const { return other <= value_; }
  
  inline binary& operator+= ( const binary& other ) { value_ += other.value_; return *this; }
  inline binary& operator+= ( const uint64_t& other ) { value_ += other; return *this; }
  inline binary& operator|= ( const binary& other ) { value_ |= other.value_; return *this; }
  inline binary& operator|= ( const uint64_t& other ) { value_ |= other; return *this; }
  uint64_t value_;
};

inline std::ostream& operator<< ( std::ostream& aStr , const binary& value ){
  std::ios_base::fmtflags f( std::cout.flags() ); 
  aStr << std::hex << std::setfill('0') << std::setw(16) << value.value_; 
  std::cout.flags( f );
  return aStr;
}

inline std::istream& operator>> ( std::istream& aStr , binary& value ){ return ( aStr >> std::hex >> value.value_ >> std::dec ); }

  
#endif