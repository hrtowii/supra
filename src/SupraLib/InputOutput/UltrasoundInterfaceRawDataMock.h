// ================================================================================================
// 
// If not explicitly stated: Copyright (C) 2016, all rights reserved,
//      Rüdiger Göbl 
//		Email r.goebl@tum.de
//      Chair for Computer Aided Medical Procedures
//      Technische Universität München
//      Boltzmannstr. 3, 85748 Garching b. München, Germany
// 
// ================================================================================================


#ifndef __ULTRASOUNDINTERFACERAWDATAMOCK_H__
#define __ULTRASOUNDINTERFACERAWDATAMOCK_H__

#ifdef HAVE_BEAMFORMER

#include <atomic>
#include <memory>
#include <mutex>

#include <AbstractInput.h>
#include <USImage.h>
#include <Container.h>

namespace supra
{
	//forward declaration
	class USRawData;

	class UltrasoundInterfaceRawDataMock : public AbstractInput<RecordObject>
	{
	public:
		UltrasoundInterfaceRawDataMock(tbb::flow::graph& graph, const std::string & nodeID);

		//Functions to be overwritten
	public:
		virtual void initializeDevice();
		virtual bool ready();

		virtual std::vector<size_t> getImageOutputPorts() { return{}; };
		virtual std::vector<size_t> getTrackingOutputPorts() { return{}; };

		virtual void freeze();
		virtual void unfreeze();
	protected:
		virtual void startAcquisition();
		//Needs to be thread safe
		virtual void configurationEntryChanged(const std::string& configKey);
		//Needs to be thread safe
		virtual void configurationChanged();

		virtual bool timerCallback();

	private:
		static constexpr size_t m_maxSequenceLength = 20;
		static constexpr size_t m_maxSequenceSizeMb = 512;

		void readConfiguration();
		void readNextFrame();

		std::string m_mockMetadataFilename;
		std::vector<std::string> m_mockDataFilenames;
		bool m_singleImage;
		bool m_streamSequenceOnce;
		int m_frequency;
		std::shared_ptr<USRawData> m_protoRawData;
		std::shared_ptr<Container<int16_t> > m_pMockData;

		std::vector<std::shared_ptr<std::ifstream>> m_mockDataStreams;
		std::vector<std::vector<char> > m_mockDataStramReadBuffers;

		std::vector<size_t> m_sequenceLengths;
		size_t m_sequenceIndex;
		size_t m_frameIndex;
		size_t m_numel;
		std::atomic_bool m_frozen;

		std::mutex m_objectMutex;
	};
}

#endif //!HAVE_BEAMFORMER

#endif //!__ULTRASOUNDINTERFACERAWDATAMOCK_H__
