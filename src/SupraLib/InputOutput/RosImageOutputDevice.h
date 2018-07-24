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

#ifndef __ROSIMAGEOUTPUTDEVICE_H___
#define __ROSIMAGEOUTPUTDEVICE_H__

#ifdef HAVE_DEVICE_ROSIMAGE_OUTPUT

#include "AbstractOutput.h"

namespace supra
{
	//forward declarations
	class USImage;
	class RosWrapper;

	class RosImageOutputDevice : public AbstractOutput
	{
	public:
		RosImageOutputDevice(tbb::flow::graph& graph, const std::string& nodeID, bool queueing);
		~RosImageOutputDevice();

		//Functions to be overwritten
	public:
		virtual void initializeOutput();
		virtual bool ready();
	protected:
		virtual void startOutput();
		//Needs to be thread safe
		virtual void stopOutput();
		//Needs to be thread safe
		virtual void configurationDone();

		virtual void writeData(std::shared_ptr<RecordObject> data);

	private:
		void addData(std::shared_ptr<const RecordObject> data);
		void addSyncRecord(std::shared_ptr<const RecordObject> _syncMessage);
		template <typename ElementType>
		void addImageTemplated(std::shared_ptr<const USImage> imageData);
		void addImage(std::shared_ptr<const RecordObject> imageData);

		bool m_isReady;

		std::string m_topic;
		std::string m_masterHost;

		std::unique_ptr<RosWrapper> m_rosWrapper;
		size_t m_publisherNoImage;
	};
}

#endif //!HAVE_DEVICE_ROSIMAGE_OUTPUT

#endif //!__ROSIMAGEOUTPUTDEVICE_H__
