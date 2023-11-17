#ifndef L1Trigger_L1THGCal_Toolkit_h
#define L1Trigger_L1THGCal_Toolkit_h

#include <iostream>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
#include <memory>



template< typename T >
inline std::ostream& print( std::ostream& aStr , const T& aArg )
{
  return aStr << std::endl;
}

template< typename T , typename X , typename... Args >
inline std::ostream& print( std::ostream& aStr , const T& aArg , const X& aMember , const Args&... aRest )
{
  aStr << std::setw(11) << aArg.*aMember << " ";
  return print( aStr , aArg , std::move( aRest )... );
}


template< typename T >
inline bool check( const T& aArg1 ,  const T& aArg2 )
{
  return true;
}

template< typename T , typename X , typename... Args >
inline bool check( const T& aArg1 ,  const T& aArg2 , const X& aMember , const Args&... aRest )
{
  if( aArg1.*aMember != aArg2.*aMember ) return false;
  return check( aArg1 ,  aArg2 , std::move( aRest )... );
}


template< typename T  > 
void check( const std::string& aTitle , const std::vector< std::unique_ptr< T > >& aLeft , const std::vector< std::unique_ptr< T > >& aRight )
{
  std::cout << aTitle << ": " << std::flush;
  // if( aLeft.size() != aRight.size() )
  // {
    // std::cout << "\nData:     " << aRight.size() << " elements     vs.     Expected: " << aLeft.size()  << " elements" << std::endl;
    // throw std::runtime_error( "Vector size mismatch" );
  // }
  
  bool lError( false );
  auto lIt1( aLeft.begin() ) , lIt2( aRight.begin() );
  for(  ; (lIt1 < aLeft.end()) && (lIt2 < aRight.end()) ; ++lIt1 , ++lIt2 )
  {
    if( not (**lIt1 == **lIt2 ) )
    {
      std::cout << "\nData:     " << **lIt2 << "Expected: " << **lIt1 << std::flush;
      lError = true;
    }
  }
  
  if( aLeft.size() != aRight.size() )
  {
    std::cout << "\nData:     " << aRight.size() << " elements     vs.     Expected: " << aLeft.size()  << " elements" << std::endl;
    
    for(  ; (lIt1 < aLeft.end()) ; ++lIt1 )  std::cout << "\nExpected: " << **lIt1 << std::flush;
    for(  ; (lIt2 < aRight.end()) ; ++lIt2 ) std::cout << "\nData:     " << **lIt2 << std::flush;
 
    
    lError = true;
  }
  
  if( lError ) throw std::runtime_error( "Vector data mismatch" );
  std::cout << "Pass" << std::endl;
}

template< typename T  > 
void check( const std::string& aTitle , const std::vector< std::shared_ptr< T > >& aLeft , const std::vector< std::shared_ptr< T > >& aRight )
{
  std::cout << aTitle << ": " << std::flush;
  // if( aLeft.size() != aRight.size() )
  // {
    // std::cout << "\nData:     " << aRight.size() << " elements     vs.     Expected: " << aLeft.size()  << " elements" << std::endl;
    // throw std::runtime_error( "Vector size mismatch" );
  // }
  
  bool lError( false );
  auto lIt1( aLeft.begin() ) , lIt2( aRight.begin() );
  for(  ; (lIt1 < aLeft.end()) && (lIt2 < aRight.end()) ; ++lIt1 , ++lIt2 )
  {
    if( not (**lIt1 == **lIt2 ) )
    {
      std::cout << "\nData:     " << **lIt2 << "Expected: " << **lIt1 << std::flush;
      lError = true;
    }
  }
  
  if( aLeft.size() != aRight.size() )
  {
    std::cout << "\nData:     " << aRight.size() << " elements     vs.     Expected: " << aLeft.size()  << " elements" << std::endl;
    
    for(  ; (lIt1 < aLeft.end()) ; ++lIt1 )  std::cout << "\nExpected: " << **lIt1 << std::flush;
    for(  ; (lIt2 < aRight.end()) ; ++lIt2 ) std::cout << "\nData:     " << **lIt2 << std::flush;
 
    
    lError = true;
  }
  
  if( lError ) throw std::runtime_error( "Vector data mismatch" );
  std::cout << "Pass" << std::endl;
}


template< typename T  > 
void columnwise_check( const std::string& aTitle , const std::vector< std::unique_ptr< T > >& aLeft , const std::vector< std::unique_ptr< T > >& aRight , const std::size_t& stride )
{
  std::cout << aTitle << ": " << std::flush;
  if( aLeft.size() != aRight.size() ) throw std::runtime_error( "Vector size mismatch" );
  
  bool lError( false );
  auto lStart1( aLeft.begin() ) , lStart2( aRight.begin() );
  for( int i(0); i!=stride; ++i, ++lStart1, ++lStart2 )
  {
    for( auto lIt1( lStart1 ) , lIt2( lStart2 ) ; lIt1 < aLeft.end() ; lIt1+=stride , lIt2+=stride )
    {
      if( not (**lIt1 == **lIt2 ) )
      {
        std::cout << "\nData:     " << **lIt2 << "Expected: " << **lIt1 << std::flush;
        lError = true;
      }
    }
    if( lError ) throw std::runtime_error( "Vector data mismatch" );
    std::cout << "Pass" << std::endl;
  }
}




template< typename T  >
inline void ReadDebugFile( std::stringstream& aLine , T& aRet )
{}

template< typename T , typename... Args >
inline void ReadDebugFile( std::stringstream& aLine , T& aRet , bool T::*aMember , const Args&... aRest )
{
  std::string lVal;
  aLine >> lVal;
  aRet.*aMember = (lVal == "TRUE");
  ReadDebugFile( aLine , aRet , std::move( aRest )... );
}

template< typename T , typename X , typename... Args >
inline void ReadDebugFile( std::stringstream& aLine , T& aRet , X& aMember , const Args&... aRest )
{
  aLine >> aRet.*aMember;
  ReadDebugFile( aLine , aRet , std::move( aRest )... ); 
}


template< typename T , typename... Args >
inline void ReadDebugFile( const std::string& aFilename , std::vector< std::unique_ptr< T > >& aRet , const Args&... aRest )
{
  std::ifstream infile( aFilename ); 
  std::string line; 
  
  for( int i(0); i != 2; ++i ) std::getline( infile, line ); // Headers + offset
  
  while ( std::getline( infile, line ) )
  { 
    std::stringstream lStr( line );
    auto lMember = std::make_unique< T >();
    ReadDebugFile( lStr , *lMember , std::move( aRest )... ); 
    aRet.push_back( move(lMember) );
  }
}

template< typename T , typename... Args >
inline void ReadDebugFile( const std::string& aFilename , std::vector< std::shared_ptr< T > >& aRet , const Args&... aRest )
{
  std::ifstream infile( aFilename ); 
  std::string line; 
  
  for( int i(0); i != 2; ++i ) std::getline( infile, line ); // Headers + offset
  
  while ( std::getline( infile, line ) )
  { 
    std::stringstream lStr( line );
    auto lMember = std::make_shared< T >();
    ReadDebugFile( lStr , *lMember , aRest... ); 
    aRet.push_back( lMember );
  }
}

template< typename T  >
inline std::vector< T >& resize( std::vector< T >&& aVector , const std::size_t& aSize )
{
  aVector.resize( aSize );
  return aVector;
}


template< typename T  >
inline std::vector< T >& resize( std::vector< T >&& aVector , const std::size_t& aOffset , const std::size_t& aSize )
{
  aVector = std::vector< T >( aVector.begin() + aOffset , aVector.begin() + aOffset + aSize );
  return aVector;
}


#endif