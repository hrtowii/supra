// ================================================================================================
// 
// If not explicitly stated: Copyright (C) 2011-2016, all rights reserved,
//      Christoph Hennersperger 
//		EmaiL christoph.hennersperger@tum.de
//      Chair for Computer Aided Medical Procedures
//      Technische Universität München
//      Boltzmannstr. 3, 85748 Garching b. München, Germany
//	and
//		Rüdiger Göbl
//		Email r.goebl@tum.de
//
// ================================================================================================

#include "USImage.h"
#include "Beamformer/USRawData.h"
#include "UltrasoundInterfaceRawDataMock.h"
#include "utilities/utility.h"

#include "Beamformer/RxBeamformerParameters.h"
#include "ContainerFactory.h"

#include <memory>

using namespace std;

namespace supra
{
	UltrasoundInterfaceRawDataMock::UltrasoundInterfaceRawDataMock(tbb::flow::graph & graph, const std::string & nodeID)
		: AbstractInput(graph, nodeID,1)
		, m_sequenceIndex(0)
		, m_frameIndex(0)
		, m_numel(0)
		, m_frozen(false)
		, m_lastFrame(false)
		, m_ready(false)
	{
		m_callFrequency.setName("RawMock");
		//Setup allowed values for parameters
		m_valueRangeDictionary.set<bool>("singleImage", { true, false }, false, "Single image");
		m_valueRangeDictionary.set<bool>("streamSequenceOnce", { true, false }, false, "Emit sequences once");
		m_valueRangeDictionary.set<double>("frequency", 0.001, 100, 5, "Frequency");
		m_valueRangeDictionary.set<string>("mockMetaDataFilename", "", "Mock meta data filename");
		m_valueRangeDictionary.set<string>("mockDataFilename", "", "Mock data filename");

		readConfiguration();
	}

	void UltrasoundInterfaceRawDataMock::initializeDevice()
	{
		if (getTimerFrequency() != m_frequency)
		{
			setUpTimer(m_frequency);
		}

		// try
		// {
		// 	m_protoRawData = RxBeamformerParameters::readMetaDataForMock(m_mockMetadataFilename);

		// 	m_numel = m_protoRawData->getNumReceivedChannels()*m_protoRawData->getNumSamples()*m_protoRawData->getNumScanlines();

		// 	// initialize m_mockDataStreams and m_sequenceLengths by getting the file sizes of all datafiles
		// 	m_mockDataStramReadBuffers.resize(m_mockDataFilenames.size());
		// 	m_mockDataStreams.resize(m_mockDataFilenames.size());
		// 	m_sequenceLengths.resize(m_mockDataFilenames.size());
		// 	for (size_t k = 0; k < m_mockDataFilenames.size(); k++)
		// 	{
		// 		// In order to maximize reading performance, the ifstream needs a large read buffer
		// 		m_mockDataStramReadBuffers[k].resize(128 * 1024, '\0');
		// 		m_mockDataStreams[k] = std::shared_ptr<std::ifstream>(new std::ifstream);
		// 		m_mockDataStreams[k]->open(m_mockDataFilenames[k], std::ifstream::ate | std::ifstream::binary);
		// 		if (!m_mockDataStreams[k]->good())
		// 		{
		// 			logging::log_error("UltrasoundInterfaceRawDataMock: Error opening mock file ", m_mockDataFilenames[k]);
		// 			throw std::runtime_error("UltrasoundInterfaceRawDataMock: Error opening mock file ");
		// 		}
		// 		m_mockDataStreams[k]->rdbuf()->pubsetbuf(m_mockDataStramReadBuffers[k].data(), m_mockDataStramReadBuffers[k].size());	
		// 		size_t filesizeBytes = m_mockDataStreams[k]->tellg();
		// 		m_mockDataStreams[k]->seekg(0);

		// 		m_sequenceLengths[k] = filesizeBytes / (m_numel * sizeof(int16_t));
		// 	}

		// 	readNextFrame();
		// 	m_ready = true;
		// }
		// catch (std::exception e)
		// {
		// 	logging::log_error("UltrasoundInterfaceRawDataMock: Caught exception preparing mock meta or mock file");
		// 	m_ready = false;
		// }
		try
{
    m_protoRawData = RxBeamformerParameters::readMetaDataForMock(m_mockMetadataFilename);
    
    // Debug the parsed values
    auto channels = m_protoRawData->getNumReceivedChannels();
    auto samples = m_protoRawData->getNumSamples();
    auto scanlines = m_protoRawData->getNumScanlines();
    
    logging::log_info("Raw parsed metadata values:");
    logging::log_info("  getNumReceivedChannels(): ", channels);
    logging::log_info("  getNumSamples(): ", samples);
    logging::log_info("  getNumScanlines(): ", scanlines);
    
    m_numel = channels * samples * scanlines;
    logging::log_info("Calculated m_numel: ", m_numel);
    
    // Sanity check against actual file size
    std::ifstream testFile(m_mockDataFilenames[0], std::ios::binary | std::ios::ate);
    if (!testFile.is_open()) {
        throw std::runtime_error("Could not open file for size check: " + m_mockDataFilenames[0]);
    }
    size_t actualFileSize = testFile.tellg();
    testFile.close();
    
    size_t expectedNumelFromFile = actualFileSize / sizeof(int16_t);
    logging::log_info("File analysis:");
    logging::log_info("  Actual file size: ", actualFileSize, " bytes (", actualFileSize/1024/1024, " MB)");
    logging::log_info("  Expected m_numel from file size: ", expectedNumelFromFile);
    
    // CRITICAL FIX: Override if parser gave unreasonable values
    if (m_numel > expectedNumelFromFile * 2 || m_numel < expectedNumelFromFile / 2) {
        logging::log_info("Parser gave unreasonable m_numel (", m_numel, "), using file-based calculation");
        m_numel = expectedNumelFromFile;
    }
    
    // Additional safety limit
    const size_t MAX_REASONABLE_NUMEL = 100000000;  // 100M elements = ~200MB
    if (m_numel > MAX_REASONABLE_NUMEL) {
        logging::log_error("m_numel still too large, capping at safety limit");
        m_numel = std::min(m_numel, expectedNumelFromFile);
    }
    
    logging::log_info("FINAL m_numel: ", m_numel, " (", (m_numel * sizeof(int16_t))/1024/1024, " MB)");
    
    // Rest of initialization...
    m_mockDataStramReadBuffers.resize(m_mockDataFilenames.size());
    m_mockDataStreams.resize(m_mockDataFilenames.size());
    m_sequenceLengths.resize(m_mockDataFilenames.size());
    
    for (size_t k = 0; k < m_mockDataFilenames.size(); k++)
    {
        m_mockDataStramReadBuffers[k].resize(128 * 1024, '\0');
        m_mockDataStreams[k] = std::shared_ptr<std::ifstream>(new std::ifstream);
        m_mockDataStreams[k]->open(m_mockDataFilenames[k], std::ifstream::ate | std::ifstream::binary);
        
        if (!m_mockDataStreams[k]->good()) {
            logging::log_error("Error opening mock file ", m_mockDataFilenames[k]);
            throw std::runtime_error("Error opening mock file");
        }
        
        m_mockDataStreams[k]->rdbuf()->pubsetbuf(m_mockDataStramReadBuffers[k].data(), m_mockDataStramReadBuffers[k].size());
        size_t filesizeBytes = m_mockDataStreams[k]->tellg();
        m_mockDataStreams[k]->seekg(0);
        m_sequenceLengths[k] = filesizeBytes / (m_numel * sizeof(int16_t));
    }
    
    readNextFrame();
    m_ready = true;
}
catch (const std::bad_alloc& e) {
    logging::log_error("BAD_ALLOC: ", e.what(), " - attempted ", (m_numel * sizeof(int16_t))/1024/1024, " MB");
    m_ready = false;
}
catch (const std::exception& e) {
    logging::log_error("Exception: ", typeid(e).name(), " - ", e.what());
    m_ready = false;
}
	}

	void UltrasoundInterfaceRawDataMock::freeze()
	{
		m_frozen = true;
	}

	void UltrasoundInterfaceRawDataMock::unfreeze()
	{
		m_frozen = false;
	}

	void UltrasoundInterfaceRawDataMock::startAcquisition()
	{
		setUpTimer(m_frequency);
		timerLoop();
	}

	void UltrasoundInterfaceRawDataMock::configurationEntryChanged(const std::string & configKey)
	{
		lock_guard<mutex> lock(m_objectMutex);
		if (configKey == "frequency")
		{
			m_frequency = m_configurationDictionary.get<double>("frequency");
			if (getTimerFrequency() != m_frequency)
			{
				setUpTimer(m_frequency);
			}
		}
		if (configKey == "singleImage")
		{
			m_singleImage = m_configurationDictionary.get<bool>("singleImage");
		}
		if (configKey == "streamSequenceOnce")
		{
			m_streamSequenceOnce = m_configurationDictionary.get<bool>("streamSequenceOnce");
		}
	}

	void UltrasoundInterfaceRawDataMock::configurationChanged()
	{
		readConfiguration();
	}

	bool UltrasoundInterfaceRawDataMock::timerCallback() {
		if (!m_frozen)
		{
			double timestamp = getCurrentTime();

			m_callFrequency.measure();
			shared_ptr<USRawData> pRawData = std::make_shared<USRawData>(
				m_protoRawData->getNumScanlines(),
				m_protoRawData->getNumElements(),
				m_protoRawData->getElementLayout(),
				m_protoRawData->getNumReceivedChannels(),
				m_protoRawData->getNumSamples(),
				m_protoRawData->getSamplingFrequency(),
				m_pMockData,
				m_protoRawData->getRxBeamformerParameters(),
				m_protoRawData->getImageProperties(),
				getCurrentTime(),
				getCurrentTime());
			addData<0>(pRawData);

			if (!m_singleImage)
			{
				if (m_lastFrame)
				{
					setRunning(false);
				}
				else
				{
					readNextFrame();
				}
			}
			m_callFrequency.measureEnd();
		}
		return getRunning();
	}

	void UltrasoundInterfaceRawDataMock::readConfiguration()
	{
		lock_guard<mutex> lock(m_objectMutex);
		//read conf values
		m_singleImage = m_configurationDictionary.get<bool>("singleImage");
		m_streamSequenceOnce = m_configurationDictionary.get<bool>("streamSequenceOnce");
		m_frequency = m_configurationDictionary.get<double>("frequency");
		m_mockMetadataFilename = m_configurationDictionary.get<string>("mockMetaDataFilename");
		m_mockDataFilenames = split(m_configurationDictionary.get<string>("mockDataFilename"), ',');
		for (auto& filename : m_mockDataFilenames)
		{
			filename = trim(filename);
		}
	}

	// void UltrasoundInterfaceRawDataMock::readNextFrame()
	// {
	// 	auto mockDataHost = make_shared<Container<int16_t> >(LocationHost, ContainerFactory::getNextStream(), m_numel);

	// 	m_mockDataStreams[m_sequenceIndex]->read(reinterpret_cast<char*>(mockDataHost->get()), m_numel * sizeof(int16_t));
	// 	m_pMockData = make_shared<Container<int16_t> >(LocationGpu, *mockDataHost);
	// 	// advance to the next image and sequence where required
	// 	m_frameIndex = (m_frameIndex + 1) % m_sequenceLengths[m_sequenceIndex];
	// 	if (m_frameIndex == 0)
	// 	{
	// 		m_mockDataStreams[m_sequenceIndex]->seekg(0);
	// 		m_sequenceIndex = (m_sequenceIndex + 1) % m_sequenceLengths.size();
	// 		if (m_sequenceIndex == 0 && m_streamSequenceOnce)
	// 		{
	// 			m_lastFrame = true;
	// 		}
	// 	}
	// }
	void UltrasoundInterfaceRawDataMock::readNextFrame()
{
    // Safety check before massive allocation
    size_t frameSizeBytes = m_numel * sizeof(int16_t);
    if (frameSizeBytes > 500 * 1024 * 1024) {  // 500 MB limit
        throw std::runtime_error("Frame size too large: " + std::to_string(frameSizeBytes/1024/1024) + " MB");
    }
    
    // Force ALL allocations to use host memory only
    auto mockDataHost = make_shared<Container<int16_t>>(LocationHost, ContainerFactory::getNextStream(), m_numel);
    m_mockDataStreams[m_sequenceIndex]->read(reinterpret_cast<char*>(mockDataHost->get()), m_numel * sizeof(int16_t));
    
	static int frameNum = 0;
    frameNum++;
    
    if (frameNum % 10 == 1) { // Save every 10th frame to avoid spam
        logging::log_info("Dumping frame ", frameNum, " to PGM");
        
        // Based on your metadata: 64 channels, 712704 samples, 5568 scanlines
        // Let's visualize one scanline (64 channels x numSamples)
        int channels = 64;
        int samples = 712704 / 5568; // samples per scanline = ~128
        int scanlineToShow = 0; // first scanline
        
        // Create a simple bitmap
        std::vector<uint8_t> imageData(channels * samples);
        
        // Find min/max for scaling
        int16_t minVal = 32767, maxVal = -32768;
        for (int i = 0; i < channels * samples; i++) {
            int16_t val = mockDataHost->get()[scanlineToShow * channels * samples + i];
            minVal = std::min(minVal, val);
            maxVal = std::max(maxVal, val);
        }
        
        logging::log_info("Data range: [", minVal, ", ", maxVal, "]");
        
        // Convert to 0-255 range
        for (int row = 0; row < samples; row++) {
            for (int col = 0; col < channels; col++) {
                int16_t val = mockDataHost->get()[scanlineToShow * channels * samples + row * channels + col];
                double normalized = (maxVal != minVal) ? 
                    static_cast<double>(val - minVal) / (maxVal - minVal) : 0.0;
                imageData[row * channels + col] = static_cast<uint8_t>(normalized * 255.0);
            }
        }
        
        // Write simple PGM file to /supra/data
        char filename[256];
        snprintf(filename, sizeof(filename), "/supra/data/frame_%06d.pgm", frameNum);
        
        FILE* f = fopen(filename, "wb");
        if (f) {
            fprintf(f, "P5\n%d %d\n255\n", channels, samples);
            fwrite(imageData.data(), 1, channels * samples, f);
            fclose(f);
            logging::log_info("Saved frame to: ", filename);
        } else {
            logging::log_error("Failed to open file: ", filename);
        }
    }
	
    // CRITICAL: Use LocationHost instead of LocationGpu
    m_pMockData = make_shared<Container<int16_t>>(LocationHost, *mockDataHost);
    
    // Advance frame logic...
    m_frameIndex = (m_frameIndex + 1) % m_sequenceLengths[m_sequenceIndex];
    if (m_frameIndex == 0) {
        m_mockDataStreams[m_sequenceIndex]->seekg(0);
        m_sequenceIndex = (m_sequenceIndex + 1) % m_sequenceLengths.size();
        if (m_sequenceIndex == 0 && m_streamSequenceOnce) {
            m_lastFrame = true;
        }
    }
	// logging::log_info("MOCK INTERFACE: Created frame data, type = ", m_pMockData ? "valid" : "NULL");
    // if (m_pMockData) {
    //     logging::log_info("MOCK INTERFACE: Frame data size = ", m_pMockData->size());
    // }
}

	bool UltrasoundInterfaceRawDataMock::ready()
	{
		return m_ready;
	}
}
