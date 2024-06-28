#include "L1Trigger/L1THGCal/interface/mappingTools/HgcConfigReader.hpp"
#include "L1Trigger/L1THGCal/interface/pugixml-src/pugixml.hpp"
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include "xercesc/util/XercesDefs.hpp"
#include <iostream>
#include <iomanip>
#include <numeric>

XERCES_CPP_NAMESPACE_USE

namespace l1thgcmapping {

  template< typename T , typename U , typename V , typename W  >
  void ReciprocalInsert( const T& a , const U& b , V& a_to_b , W& b_to_a )
  {
    a_to_b.insert( {{ a , b }} );
    b_to_a.insert( {{ b , a }} );
  }

  template< typename T , typename U , typename V , typename W  >
  void ReciprocalInsert_XercesC( const DOMDocument* doc , const std::string& aXpath , const T& aChildIdFn , const U& aParentIdFn , V& a_to_b , W& b_to_a )
  {

    // Loop nodes
    DOMNodeList* nodeList = doc->getElementsByTagName(XMLString::transcode(aXpath.c_str()));
    const XMLSize_t length = nodeList->getLength();
    for (XMLSize_t i = 0; i < length; ++i) {
      DOMNode* node = nodeList->item(i);

      // Get node
      DOMElement* element = dynamic_cast<DOMElement*>(node);

      // Get parent
      DOMNode* parent = node->getParentNode();
      DOMElement* element_parent = dynamic_cast<DOMElement*>(parent);

      // Insert into map
      ReciprocalInsert( aChildIdFn(element) , aParentIdFn(element_parent) , a_to_b , b_to_a );
    }
  }

  const auto getNode = [] ( const DOMElement* element, const std::string& nodeName ) { 
    const XMLCh* idAttr = XMLString::transcode(nodeName.c_str());
    const XMLCh* idValue = element->getAttribute(idAttr);
    char* idStr = XMLString::transcode(idValue);
    return idStr;
  };

  const auto DecValue_XercesC = []( const char* value ){
    return std::stoul(value, nullptr , 10 );
  };
  const auto HexValue_XercesC = []( const char* value ){
    return std::stoul(value, nullptr , 16 );
  };

  const auto DecId_XercesC = []( const DOMElement* element ){
    return std::stoul(getNode(element, "id"), nullptr , 10 );
  };
  const auto HexId_XercesC = []( const DOMElement* element ){
    return std::stoul(getNode(element, "id"), nullptr , 16 );
  };
  const auto HexHref_XercesC = []( const DOMElement* element ){
    return std::stoul(getNode(element, "href"), nullptr , 16 );
  };

  DOMDocument* XmlOpen_XercesC( const std::string& aFilename )
  {
      try
      {
          XMLPlatformUtils::Initialize();

          XercesDOMParser parser;
          parser.setValidationScheme(XercesDOMParser::Val_Auto);
          parser.setDoNamespaces(true);

          parser.parse(aFilename.c_str());

          if (parser.getErrorCount() != 0)
          {
              throw std::runtime_error("Unable to parse XML file '" + aFilename + "'");
          }

          return parser.adoptDocument();
      }
      catch (const XMLException& e)
      {
          char* message = XMLString::transcode(e.getMessage());
          std::string errorMsg(message);
          XMLString::release(&message);

          throw std::runtime_error("Error initializing Xerces-C: " + errorMsg);
      }
  }

  std::string getChildXMLFileName( const DOMDocument* doc )
  {
    const auto srcAttr = doc->getDocumentElement()->getAttribute(XMLString::transcode("Src"));
    char* srcValue = XMLString::transcode(srcAttr);
    return std::string(srcValue);
  }


  void OpenGeometry( const std::string& aFilename , Maps& aMaps )
  {
    
    auto doc_xerces = XmlOpen_XercesC( aFilename.c_str() );
    ReciprocalInsert_XercesC( doc_xerces, "Module" , HexId_XercesC , HexId_XercesC , aMaps.module_to_motherboard , aMaps.motherboard_to_module );

    DOMNodeList* mbList = doc_xerces->getElementsByTagName(XMLString::transcode("Motherboard"));
    const XMLSize_t length = mbList->getLength();
    unsigned int totalNLPGBTs = 0;
    unsigned int iLPGBT = 0;
    for (XMLSize_t i = 0; i < length; ++i) {
      // Get motherboard ID and number of trigger lpgbts
      DOMNode* mb = mbList->item(i);
      DOMElement* element = dynamic_cast<DOMElement*>(mb);
      const auto mbID = HexId_XercesC(element);
      const unsigned int nTriggerLPGBTs = std::stoul(getNode(element, "TriggerLpGbts"), nullptr , 10 );

      totalNLPGBTs += nTriggerLPGBTs;
      std::vector<Maps::tID> lpgbts(nTriggerLPGBTs);
      std::iota(lpgbts.begin(), lpgbts.end(), iLPGBT);
      aMaps.motherboard_to_trigLPGBTs.insert( {{ mbID, lpgbts }} );
      iLPGBT += nTriggerLPGBTs;
      aMaps.motherboard_to_nTrigLPGBT.insert( {{ mbID, nTriggerLPGBTs }} );
    }
    std::cout << "Total N LPGBTS : " << totalNLPGBTs << std::endl;
  }

  void OpenRegions( const std::string& aFilename , Maps& aMaps )
  {   
    auto doc_xerces = XmlOpen_XercesC( aFilename.c_str() );
    OpenGeometry( getChildXMLFileName(doc_xerces), aMaps );

    ReciprocalInsert_XercesC( doc_xerces, "Motherboard" , HexHref_XercesC , HexId_XercesC , aMaps.motherboard_to_region , aMaps.region_to_motherboard );
  }

  void OpenS1( const std::string& aFilename , Maps& aMaps )
  {
    auto doc_xerces = XmlOpen_XercesC( aFilename.c_str() );
    OpenRegions( getChildXMLFileName(doc_xerces), aMaps );
    ReciprocalInsert_XercesC( doc_xerces, "Region" , HexHref_XercesC , HexId_XercesC , aMaps.region_to_stage1 , aMaps.stage1_to_region );
    // std::cout << "Number of entries in region_to_stage1 map after xerces-C : " << aMaps.region_to_stage1.size() << " " << aMaps.region_to_stage1.begin()->first << " " << aMaps.region_to_stage1.begin()->second << std::endl;
  }

  void OpenChannelAllocation( const std::string& aFilename , Maps& aMaps )
  {
    auto doc_xerces = XmlOpen_XercesC( aFilename.c_str() );
    OpenS1( getChildXMLFileName(doc_xerces), aMaps );

    // Loop channels
    DOMNodeList* channelList = doc_xerces->getElementsByTagName(XMLString::transcode("Channel"));
    const XMLSize_t nChannels = channelList->getLength();
    for (XMLSize_t iChannel = 0; iChannel < nChannels; ++iChannel) {

      DOMNode* channel = channelList->item(iChannel);
      auto ChannelId_XercesC = HexId_XercesC(dynamic_cast<DOMElement*>(channel));

      // Loop frames for this channel
      DOMNodeList* frameList = channel->getChildNodes();
      const XMLSize_t nFrames = frameList->getLength();
      for (XMLSize_t iFrame = 0; iFrame < nFrames; ++iFrame) {

        DOMNode* node = frameList->item(iFrame);
        if ( node->getNodeType() != DOMNode::ELEMENT_NODE ) continue;

        DOMElement* frame = dynamic_cast<DOMElement*>(frameList->item(iFrame));
        const XMLCh* idAttr = XMLString::transcode("Module");

        uint32_t moduleid = 0;
        if ( frame->hasAttribute( XMLString::transcode("Module") ) ) moduleid = HexValue_XercesC(getNode(frame,"Module"));
        else if ( frame->hasAttribute( XMLString::transcode("Motherboard") ) ) moduleid = HexValue_XercesC(getNode(frame,"Motherboard"));
        else continue;

        auto s1ChannelFrame = std::make_tuple( ChannelId_XercesC , DecId_XercesC( frame ) );
        auto moduleColIndex = std::make_tuple( moduleid , DecValue_XercesC( getNode(frame,"column") ) , DecValue_XercesC( getNode(frame,"index") ) );
        ReciprocalInsert( s1ChannelFrame , moduleColIndex , aMaps.s1ChannelFrame_to_moduleColIndex , aMaps.moduleColIndex_to_s1ChannelFrame );  
      }
    }
  }

  void OpenBackendMapping( const std::string& aFilename , Maps& aMaps )
  {
    auto doc_xerces = XmlOpen_XercesC( aFilename.c_str() );
    OpenChannelAllocation( getChildXMLFileName(doc_xerces) , aMaps ); 

    DOMNodeList* nodeList = doc_xerces->getElementsByTagName(XMLString::transcode("fibre"));
    const XMLSize_t length = nodeList->getLength();
    for (XMLSize_t i = 0; i < length; ++i) {
      // Get node
      DOMNode* node = nodeList->item(i);

      DOMElement* element = dynamic_cast<DOMElement*>(node);
      auto belink  = HexValue_XercesC( getNode(element, "id") );
      auto reflink = HexValue_XercesC( getNode(element, "ref-fibre") );
      auto s1out   = HexValue_XercesC( getNode(element, "s1-output") );
      auto s2in    = HexValue_XercesC( getNode(element, "s2-input") );

      aMaps.belink_to_reflink.insert( {{ belink , reflink }} );
      ReciprocalInsert(  s1out , belink , aMaps.s1out_to_belink , aMaps.belink_to_s1out );
      ReciprocalInsert(  belink , s2in , aMaps.belink_to_s2in , aMaps.s2in_to_belink );   
    }
  }
}
