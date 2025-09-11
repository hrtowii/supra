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

#include "previewBuilderQt.h"
#include "QImagePreviewReciever.h"
#include "QTrackerPreviewReciever.h"
#include "Beamformer/USRawData.h"
#include <QDir>
#include <QString>
#include <cfloat>  // for FLT_MAX
#include <algorithm>
#ifdef HAVE_CAMPVIS
#include "CampvisPreviewReciever.h"
#endif

namespace supra
{
	using namespace ::tbb::flow;
	using namespace ::std;

	previewBuilderQT::previewBuilderQT(graph & g, const std::string & nodeID, string name, QWidget* parentWidget, QVBoxLayout* targetLayout, QSize imageMaxSize, bool linearInterpolation)
		: QObject()
		, AbstractNode(nodeID, false)
		, m_nodeIn(g, 1,
			[this](const shared_ptr<RecordObject> & inMessage)
	{ processRecordObject(inMessage); })
		, m_imageMaxSize(imageMaxSize)
		, m_linearInterpolation(linearInterpolation)
		, m_haveImagePreview(false)
		, m_haveTrackingPreview(false)
		, m_parentWidget(parentWidget)
		, m_targetLayout(targetLayout)
		, m_name(name)
		, m_layerToShow(0)
	{
#ifdef HAVE_CAMPVIS
		m_pCampvisPreview = nullptr;
#endif
		m_pQimagePreview = nullptr;
		m_pTrackerPreview = nullptr;

		addTrackingPreviewWidget();
		addImagePreviewWidget();
	}

	void previewBuilderQT::removePreviews()
	{
		if (m_pTrackerPreview)
		{
			m_pTrackerPreview->setParent(nullptr);
			delete m_pTrackerPreview;
			m_pTrackerPreview = nullptr;
		}
		if (m_pQimagePreview)
		{
			m_pQimagePreview->setParent(nullptr);
			delete m_pQimagePreview;
			m_pQimagePreview = nullptr;
		}
#ifdef HAVE_CAMPVIS
		if (m_pCampvisPreview)
		{
			m_pCampvisPreview->setParent(nullptr);
			delete m_pCampvisPreview;
			m_pCampvisPreview = nullptr;
		}
#endif
	}

	void previewBuilderQT::setPreviewSize(QSize imageMaxSize)
	{
		m_imageMaxSize = imageMaxSize;
	}

	void previewBuilderQT::setLinearInterpolation(bool linearInterpolation)
	{
		m_linearInterpolation = linearInterpolation;
	}

// 	void previewBuilderQT::processRecordObject(const shared_ptr<RecordObject> inMessage)
// 	{
// 		if (inMessage)
// 		{
// 			if (inMessage->getType() == TypeUSImage)
// 			{
// 				auto inImage = std::dynamic_pointer_cast<USImage>(inMessage);
// 				logging::log_error_if(!inImage, "Error casting a record object to USImage, although the type was 'TypeUSImage'");
// 				if (inImage)
// 				{
// 					bool is2D = inImage->getSize().z == 1;
					
// 					bool useCampVis = !is2D;
// #ifdef HAVE_CAMPVIS
// 					if (useCampVis)
// 					{
// 						emit previewReadyObject(inMessage);
// 					}
// #else
// 					useCampVis = false;
// #endif
// 					if (!useCampVis)
// 					{
// 						shared_ptr<QImage> qtimage;
// 						tuple<double, double, bool> worldSize;
// 						m_layerToShow = m_layerToShow % inImage->getSize().z;

// 						qtimage = std::make_shared<QImage>(
// 							static_cast<int>(inImage->getSize().x),
// 							static_cast<int>(inImage->getSize().y),
// 							QImage::Format_Grayscale8);

// 						if (inImage->getDataType() == TypeUint8)
// 						{
// 							auto inImageData = inImage->getData<uint8_t>();
// 							if (!inImageData->isHost())
// 							{
// 								inImageData = make_shared<Container<uint8_t> >(LocationHost, *inImageData);
// 							}
// 							for (size_t row = 0; row < inImage->getSize().y; row++)
// 							{
// 								std::memcpy(qtimage->scanLine(static_cast<int>(row)),
// 									inImageData->get() + row*inImage->getSize().x + m_layerToShow * inImage->getSize().x*inImage->getSize().y,
// 									inImage->getSize().x * sizeof(uint8_t));
// 							}
// 						}
// 						else if (inImage->getDataType() == TypeInt16)
// 						{
// 							auto inImageData = inImage->getData<int16_t>();
// 							if (!inImageData->isHost())
// 							{
// 								inImageData = make_shared<Container<int16_t> >(LocationHost, *inImageData);
// 							}
// 							for (size_t row = 0; row < inImage->getSize().y; row++)
// 							{
// 								uchar* destRow = qtimage->scanLine(static_cast<int>(row));
// 								const int16_t* srcRow = inImageData->get() + row*inImage->getSize().x + m_layerToShow * inImage->getSize().x*inImage->getSize().y;
// 								for (size_t col = 0; col < inImage->getSize().x; col++)
// 								{
// 									destRow[col] = static_cast<uint8_t>(min(static_cast<double>(abs(srcRow[col])), 255.0));
// 								}
// 							}
// 						}
// 						else if (inImage->getDataType() == TypeFloat)
// 						{
// 							auto inImageData = inImage->getData<float>();
// 							if (!inImageData->isHost())
// 							{
// 								inImageData = make_shared<Container<float> >(LocationHost, *inImageData);
// 							}
// 							for (size_t row = 0; row < inImage->getSize().y; row++)
// 							{
// 								uchar* destRow = qtimage->scanLine(static_cast<int>(row));
// 								const float* srcRow = inImageData->get() + row*inImage->getSize().x + m_layerToShow * inImage->getSize().x*inImage->getSize().y;
// 								for (size_t col = 0; col < inImage->getSize().x; col++)
// 								{
// 									destRow[col] = static_cast<uint8_t>(min(static_cast<double>(abs(srcRow[col])), 255.0));
// 								}
// 							}
// 						}
// 						worldSize = computeWorldSize(inImage);
// 						m_layerToShow++;

// 						double worldWidth = get<0>(worldSize);
// 						double worldHeight = get<1>(worldSize);
// 						bool keepRatio = get<2>(worldSize);

// 						int imageWidth;
// 						int imageHeight;

// 						if (keepRatio)
// 						{
// 							if ((worldWidth / worldHeight) > (static_cast<double>(m_imageMaxSize.width()) / m_imageMaxSize.height()))
// 							{
// 								imageWidth = m_imageMaxSize.width();
// 								imageHeight = m_imageMaxSize.width() / worldWidth*worldHeight;
// 							}
// 							else {
// 								imageWidth = m_imageMaxSize.height() / worldHeight*worldWidth;
// 								imageHeight = m_imageMaxSize.height();
// 							}
// 						}
// 						else {
// 							imageWidth = m_imageMaxSize.width();
// 							imageHeight = m_imageMaxSize.height();
// 						}

// 						*qtimage = qtimage->scaled(
// 							imageWidth,
// 							imageHeight,
// 							Qt::IgnoreAspectRatio,
// 							(m_linearInterpolation ? Qt::SmoothTransformation : Qt::FastTransformation));

// 						emit previewReadyImage(qtimage);
// 					}
// 				}
// 			}
// 			else if (inMessage->getType() == TypeUSRawData)
// 			{
// 				auto inRawData = std::dynamic_pointer_cast<USRawData>(inMessage);
// 				logging::log_error_if(!inRawData, "Error casting a record object to USRawData, although the type was 'TypeUSRawData'");
// 				if (inRawData)
// 				{
// 					shared_ptr<QImage> qtimage;
// 					tuple<double, double, bool> worldSize;
// 					m_layerToShow = m_layerToShow % inRawData->getNumScanlines();

// 					auto numChannels = inRawData->getNumReceivedChannels();
// 					auto numSamples = inRawData->getNumSamples();

// 					qtimage = std::make_shared<QImage>(
// 						static_cast<int>(numChannels),
// 						static_cast<int>(numSamples),
// 						QImage::Format_Grayscale8);

// 					if (inRawData->getDataType() == TypeUint8)
// 					{
// 						auto inImageData = inRawData->getData<uint8_t>();
// 						if (!inImageData->isHost())
// 						{
// 							inImageData = make_shared<Container<uint8_t> >(LocationHost, *inImageData);
// 						}
// 						for (size_t row = 0; row < numSamples; row++)
// 						{
// 							uchar* destRow = qtimage->scanLine(static_cast<int>(row));
// 							for (size_t col = 0; col < numChannels; col++)
// 							{
// 								auto val = inImageData->get()[row + col*numSamples + m_layerToShow * numChannels*numSamples];
// 								destRow[col] = static_cast<uint8_t>(min(static_cast<double>(abs(val)), 255.0));
// 							}
// 						}
// 					}
// 					else if (inRawData->getDataType() == TypeInt16)
// 					{
// 						auto inImageData = inRawData->getData<int16_t>();
// 						if (!inImageData->isHost())
// 						{
// 							inImageData = make_shared<Container<int16_t> >(LocationHost, *inImageData);
// 						}
// 						for (size_t row = 0; row < numSamples; row++)
// 						{
// 							uchar* destRow = qtimage->scanLine(static_cast<int>(row));
// 							for (size_t col = 0; col < numChannels; col++)
// 							{
// 								auto val = inImageData->get()[row + col*numSamples + m_layerToShow * numChannels*numSamples];
// 								destRow[col] = static_cast<uint8_t>(min(static_cast<double>(abs(val)), 255.0));
// 							}
// 						}
// 					}
// 					else if (inRawData->getDataType() == TypeFloat)
// 					{
// 						auto inImageData = inRawData->getData<float>();
// 						if (!inImageData->isHost())
// 						{
// 							inImageData = make_shared<Container<float> >(LocationHost, *inImageData);
// 						}
// 						for (size_t row = 0; row < numSamples; row++)
// 						{
// 							uchar* destRow = qtimage->scanLine(static_cast<int>(row));
// 							for (size_t col = 0; col < numChannels; col++)
// 							{
// 								auto val = inImageData->get()[row + col*numSamples + m_layerToShow * numChannels*numSamples];
// 								destRow[col] = static_cast<uint8_t>(min(static_cast<double>(abs(val)), 255.0));
// 							}
// 						}
// 					}
// 					m_layerToShow++;

// 					int imageWidth = static_cast<int>(m_imageMaxSize.width()/384.0*numChannels);
// 					int imageHeight = m_imageMaxSize.height();

// 					*qtimage = qtimage->scaled(
// 						imageWidth,
// 						imageHeight,
// 						Qt::IgnoreAspectRatio,
// 						(m_linearInterpolation ? Qt::SmoothTransformation : Qt::FastTransformation));

// 					emit previewReadyImage(qtimage);
// 				}
// 			}
// 			else if (inMessage->getType() == TypeTrackerDataSet)
// 			{
// 				auto inTrackerData = std::dynamic_pointer_cast<TrackerDataSet>(inMessage);
// 				logging::log_error_if(!inTrackerData, "Error casting a record object to TrackerDataSet, although the type was 'TypeTrackerDataSet'");
// 				if (inTrackerData)
// 				{
// 					emit previewReadyTracking(inTrackerData);
// 				}
// 			}
// 		}
// 	}
void previewBuilderQT::processRecordObject(const shared_ptr<RecordObject> inMessage)
{
	// Create output directory if it doesn't exist
	static bool directoryCreated = false;
	if (!directoryCreated) {
		QDir dir;
		if (!dir.exists("/supra/images")) {
			bool success = dir.mkpath("/supra/images");
			logging::log_info("Created output directory /supra/images: ", success ? "SUCCESS" : "FAILED");
		}
		directoryCreated = true;
	}

	static int frameCounter = 0;
	frameCounter++;

	if (!inMessage) {
		logging::log_info("processRecordObject: Received null message");
		return;
	}

	logging::log_info("=== Processing Record Object ===");
	logging::log_info("Frame counter: ", frameCounter);
	logging::log_info("Message type: ", inMessage->getType());
	logging::log_info("Message timestamp: ", inMessage->getSyncTimestamp());

	if (inMessage->getType() == TypeUSImage)
	{
		logging::log_info("Processing USImage data");
		auto inImage = std::dynamic_pointer_cast<USImage>(inMessage);
		logging::log_error_if(!inImage, "Error casting a record object to USImage, although the type was 'TypeUSImage'");
		
		if (inImage)
		{
			// Log image properties
			auto imageSize = inImage->getSize();
			logging::log_info("USImage properties:");
			logging::log_info("  Size: ", imageSize.x, "x", imageSize.y, "x", imageSize.z);
			logging::log_info("  Data type: ", inImage->getDataType());
			logging::log_info("  Total elements: ", imageSize.x * imageSize.y * imageSize.z);

			bool is2D = imageSize.z == 1;
			logging::log_info("  Is 2D image: ", is2D ? "YES" : "NO");
			
			bool useCampVis = !is2D;
#ifdef HAVE_CAMPVIS
			if (useCampVis) {
				logging::log_info("Using CampVis for 3D visualization");
				emit previewReadyObject(inMessage);
			}
#else
			useCampVis = false;
			logging::log_info("CampVis not available, using Qt visualization");
#endif
			
			if (!useCampVis)
			{
				shared_ptr<QImage> qtimage;
				tuple<double, double, bool> worldSize;
				m_layerToShow = m_layerToShow % imageSize.z;

				logging::log_info("Creating QImage for layer: ", m_layerToShow);
				qtimage = std::make_shared<QImage>(
					static_cast<int>(imageSize.x),
					static_cast<int>(imageSize.y),
					QImage::Format_Grayscale8);

				if (inImage->getDataType() == TypeUint8)
				{
					logging::log_info("Processing TypeUint8 data");
					auto inImageData = inImage->getData<uint8_t>();
					
					if (!inImageData->isHost()) {
						logging::log_info("Converting GPU data to host memory");
						inImageData = make_shared<Container<uint8_t> >(LocationHost, *inImageData);
					}
					
					logging::log_info("Copying uint8 data to QImage");
					for (size_t row = 0; row < imageSize.y; row++) {
						std::memcpy(qtimage->scanLine(static_cast<int>(row)),
							inImageData->get() + row * imageSize.x + m_layerToShow * imageSize.x * imageSize.y,
							imageSize.x * sizeof(uint8_t));
					}
				}
				else if (inImage->getDataType() == TypeInt16)
				{
					logging::log_info("Processing TypeInt16 data");
					auto inImageData = inImage->getData<int16_t>();
					
					if (!inImageData->isHost()) {
						logging::log_info("Converting GPU data to host memory");
						inImageData = make_shared<Container<int16_t> >(LocationHost, *inImageData);
					}
					
					logging::log_info("Converting int16 data to uint8 for QImage");
					int16_t minVal = 32767, maxVal = -32768;
					
					// Find min/max for better contrast
					for (size_t i = 0; i < imageSize.x * imageSize.y; i++) {
						const int16_t* srcData = inImageData->get() + m_layerToShow * imageSize.x * imageSize.y;
						minVal = std::min(minVal, srcData[i]);
						maxVal = std::max(maxVal, srcData[i]);
					}
					
					logging::log_info("  Data range: [", minVal, ", ", maxVal, "]");
					
					for (size_t row = 0; row < imageSize.y; row++) {
						uchar* destRow = qtimage->scanLine(static_cast<int>(row));
						const int16_t* srcRow = inImageData->get() + row * imageSize.x + m_layerToShow * imageSize.x * imageSize.y;
						for (size_t col = 0; col < imageSize.x; col++) {
							// Improved scaling for better contrast
							double normalized = (maxVal != minVal) ? 
								static_cast<double>(abs(srcRow[col]) - abs(minVal)) / (abs(maxVal) - abs(minVal)) : 0.0;
							destRow[col] = static_cast<uint8_t>(std::min(normalized * 255.0, 255.0));
						}
					}
				}
				else if (inImage->getDataType() == TypeFloat)
				{
					logging::log_info("Processing TypeFloat data");
					auto inImageData = inImage->getData<float>();
					
					if (!inImageData->isHost()) {
						logging::log_info("Converting GPU data to host memory");
						inImageData = make_shared<Container<float> >(LocationHost, *inImageData);
					}
					
					logging::log_info("Converting float data to uint8 for QImage");
					float minVal = FLT_MAX, maxVal = -FLT_MAX;
					
					// Find min/max for better contrast
					for (size_t i = 0; i < imageSize.x * imageSize.y; i++) {
						const float* srcData = inImageData->get() + m_layerToShow * imageSize.x * imageSize.y;
						minVal = std::min(minVal, srcData[i]);
						maxVal = std::max(maxVal, srcData[i]);
					}
					
					logging::log_info("  Data range: [", minVal, ", ", maxVal, "]");
					
					for (size_t row = 0; row < imageSize.y; row++) {
						uchar* destRow = qtimage->scanLine(static_cast<int>(row));
						const float* srcRow = inImageData->get() + row * imageSize.x + m_layerToShow * imageSize.x * imageSize.y;
						for (size_t col = 0; col < imageSize.x; col++) {
							// Improved scaling for better contrast
							double normalized = (maxVal != minVal) ? 
								static_cast<double>(abs(srcRow[col]) - abs(minVal)) / (abs(maxVal) - abs(minVal)) : 0.0;
							destRow[col] = static_cast<uint8_t>(std::min(normalized * 255.0, 255.0));
						}
					}
				}
				else {
					logging::log_error("Unsupported data type for USImage: ", inImage->getDataType());
				}

				// Save image to file
				QString filename = QString("/supra/images/usimage_frame_%1_layer_%2.png")
					.arg(frameCounter, 6, 10, QChar('0'))
					.arg(m_layerToShow, 3, 10, QChar('0'));
				
				bool saveSuccess = qtimage->save(filename, "PNG");
				logging::log_info("Saved USImage to file: ", filename.toStdString(), " - ", saveSuccess ? "SUCCESS" : "FAILED");

				worldSize = computeWorldSize(inImage);
				m_layerToShow++;

				double worldWidth = get<0>(worldSize);
				double worldHeight = get<1>(worldSize);
				bool keepRatio = get<2>(worldSize);

				logging::log_info("World size: ", worldWidth, "x", worldHeight, ", keep ratio: ", keepRatio);

				int imageWidth;
				int imageHeight;

				if (keepRatio) {
					if ((worldWidth / worldHeight) > (static_cast<double>(m_imageMaxSize.width()) / m_imageMaxSize.height())) {
						imageWidth = m_imageMaxSize.width();
						imageHeight = static_cast<int>(m_imageMaxSize.width() / worldWidth * worldHeight);
					} else {
						imageWidth = static_cast<int>(m_imageMaxSize.height() / worldHeight * worldWidth);
						imageHeight = m_imageMaxSize.height();
					}
				} else {
					imageWidth = m_imageMaxSize.width();
					imageHeight = m_imageMaxSize.height();
				}

				logging::log_info("Scaling image to: ", imageWidth, "x", imageHeight);

				*qtimage = qtimage->scaled(
					imageWidth,
					imageHeight,
					Qt::IgnoreAspectRatio,
					(m_linearInterpolation ? Qt::SmoothTransformation : Qt::FastTransformation));

				logging::log_info("Emitting ReadyImage signal");
				emit previewReadyImage(qtimage);
			}
		}
	}
	else if (inMessage->getType() == TypeUSRawData)
	{
		logging::log_info("Processing USRawData");
		auto inRawData = std::dynamic_pointer_cast<USRawData>(inMessage);
		logging::log_error_if(!inRawData, "Error casting a record object to USRawData, although the type was 'TypeUSRawData'");
		
		if (inRawData)
		{
			// Log raw data properties
			auto numChannels = inRawData->getNumReceivedChannels();
			auto numSamples = inRawData->getNumSamples();
			auto numScanlines = inRawData->getNumScanlines();
			
			logging::log_info("USRawData properties:");
			logging::log_info("  Channels: ", numChannels);
			logging::log_info("  Samples: ", numSamples);
			logging::log_info("  Scanlines: ", numScanlines);
			logging::log_info("  Data type: ", inRawData->getDataType());
			logging::log_info("  Total elements: ", numChannels * numSamples * numScanlines);

			shared_ptr<QImage> qtimage;
			m_layerToShow = m_layerToShow % numScanlines;

			logging::log_info("Creating QImage for scanline: ", m_layerToShow);
			qtimage = std::make_shared<QImage>(
				static_cast<int>(numChannels),
				static_cast<int>(numSamples),
				QImage::Format_Grayscale8);

			if (inRawData->getDataType() == TypeUint8)
			{
				logging::log_info("Processing TypeUint8 raw data");
				auto inImageData = inRawData->getData<uint8_t>();
				
				if (!inImageData->isHost()) {
					logging::log_info("Converting GPU data to host memory");
					inImageData = make_shared<Container<uint8_t> >(LocationHost, *inImageData);
				}
				
				for (size_t row = 0; row < numSamples; row++) {
					uchar* destRow = qtimage->scanLine(static_cast<int>(row));
					for (size_t col = 0; col < numChannels; col++) {
						auto val = inImageData->get()[row + col * numSamples + m_layerToShow * numChannels * numSamples];
						destRow[col] = static_cast<uint8_t>(std::min(static_cast<double>(abs(val)), 255.0));
					}
				}
			}
			else if (inRawData->getDataType() == TypeInt16)
			{
				logging::log_info("Processing TypeInt16 raw data");
				auto inImageData = inRawData->getData<int16_t>();
				
				if (!inImageData->isHost()) {
					logging::log_info("Converting GPU data to host memory");
					inImageData = make_shared<Container<int16_t> >(LocationHost, *inImageData);
				}
				
				// Find min/max for better visualization
				int16_t minVal = 32767, maxVal = -32768;
				for (size_t row = 0; row < numSamples; row++) {
					for (size_t col = 0; col < numChannels; col++) {
						auto val = inImageData->get()[row + col * numSamples + m_layerToShow * numChannels * numSamples];
						minVal = std::min(minVal, val);
						maxVal = std::max(maxVal, val);
					}
				}
				
				logging::log_info("  Raw data range: [", minVal, ", ", maxVal, "]");
				
				for (size_t row = 0; row < numSamples; row++) {
					uchar* destRow = qtimage->scanLine(static_cast<int>(row));
					for (size_t col = 0; col < numChannels; col++) {
						auto val = inImageData->get()[row + col * numSamples + m_layerToShow * numChannels * numSamples];
						// Improved scaling
						double normalized = (maxVal != minVal) ? 
							static_cast<double>(abs(val) - abs(minVal)) / (abs(maxVal) - abs(minVal)) : 0.0;
						destRow[col] = static_cast<uint8_t>(std::min(normalized * 255.0, 255.0));
					}
				}
			}
			else if (inRawData->getDataType() == TypeFloat)
			{
				logging::log_info("Processing TypeFloat raw data");
				auto inImageData = inRawData->getData<float>();
				
				if (!inImageData->isHost()) {
					logging::log_info("Converting GPU data to host memory");
					inImageData = make_shared<Container<float> >(LocationHost, *inImageData);
				}
				
				// Find min/max for better visualization
				float minVal = FLT_MAX, maxVal = -FLT_MAX;
				for (size_t row = 0; row < numSamples; row++) {
					for (size_t col = 0; col < numChannels; col++) {
						auto val = inImageData->get()[row + col * numSamples + m_layerToShow * numChannels * numSamples];
						minVal = std::min(minVal, val);
						maxVal = std::max(maxVal, val);
					}
				}
				
				logging::log_info("  Raw data range: [", minVal, ", ", maxVal, "]");
				
				for (size_t row = 0; row < numSamples; row++) {
					uchar* destRow = qtimage->scanLine(static_cast<int>(row));
					for (size_t col = 0; col < numChannels; col++) {
						auto val = inImageData->get()[row + col * numSamples + m_layerToShow * numChannels * numSamples];
						// Improved scaling
						double normalized = (maxVal != minVal) ? 
							static_cast<double>(abs(val) - abs(minVal)) / (abs(maxVal) - abs(minVal)) : 0.0;
						destRow[col] = static_cast<uint8_t>(std::min(normalized * 255.0, 255.0));
					}
				}
			}
			else {
				logging::log_error("Unsupported data type for USRawData: ", inRawData->getDataType());
			}

			// Save raw data image to file
			QString filename = QString("/supra/images/rawdata_frame_%1_scanline_%2.png")
				.arg(frameCounter, 6, 10, QChar('0'))
				.arg(m_layerToShow, 3, 10, QChar('0'));
			
			bool saveSuccess = qtimage->save(filename, "PNG");
			logging::log_info("Saved USRawData to file: ", filename.toStdString(), " - ", saveSuccess ? "SUCCESS" : "FAILED");

			m_layerToShow++;

			int imageWidth = static_cast<int>(m_imageMaxSize.width() / 384.0 * numChannels);
			int imageHeight = m_imageMaxSize.height();

			logging::log_info("Scaling raw data image to: ", imageWidth, "x", imageHeight);

			*qtimage = qtimage->scaled(
				imageWidth,
				imageHeight,
				Qt::IgnoreAspectRatio,
				(m_linearInterpolation ? Qt::SmoothTransformation : Qt::FastTransformation));

			logging::log_info("Emitting previewReadyImage signal for raw data");
			emit previewReadyImage(qtimage);
		}
	}
	else if (inMessage->getType() == TypeTrackerDataSet)
	{
		logging::log_info("Processing TrackerDataSet");
		auto inTrackerData = std::dynamic_pointer_cast<TrackerDataSet>(inMessage);
		logging::log_error_if(!inTrackerData, "Error casting a record object to TrackerDataSet, although the type was 'TypeTrackerDataSet'");
		
		if (inTrackerData) {
			// logging::log_info("  Number of tracking entries: ", inTrackerData->getSize());
			logging::log_info("Emitting previewReadyTracking signal");
			emit previewReadyTracking(inTrackerData);
		}
	}
	else 
	{
		logging::log_info("Unknown message type received: ", inMessage->getType());
		logging::log_info("Known types are: TypeUSImage=", TypeUSImage, ", TypeUSRawData=", TypeUSRawData, ", TypeTrackerDataSet=", TypeTrackerDataSet);
	}

	logging::log_info("=== Finished Processing Record Object ===");
}
	
void previewBuilderQT::addImagePreviewWidget()
	{
		if (!m_haveImagePreview)
		{
#ifdef HAVE_CAMPVIS
			m_pCampvisPreview = new CampvisPreviewReciever(m_parentWidget, m_targetLayout, m_name);
			qRegisterMetaType<shared_ptr<RecordObject>>("std::shared_ptr<RecordObject>");
			connect(this, SIGNAL(previewReadyObject(const std::shared_ptr<RecordObject>)), m_pCampvisPreview, SLOT(previewReadyImage(const std::shared_ptr<RecordObject>)));
#endif
			//create the QImagePreviewReciever
			m_pQimagePreview = new QImagePreviewReciever(m_parentWidget, m_targetLayout, QString::fromStdString(m_name));
			qRegisterMetaType<shared_ptr<QImage>>("std::shared_ptr<QImage>");
			connect(this, SIGNAL(previewReadyImage(const std::shared_ptr<QImage>)), m_pQimagePreview, SLOT(previewReadyImage(const std::shared_ptr<QImage>)));

			m_haveImagePreview = true;
		}
	}

	void previewBuilderQT::addTrackingPreviewWidget()
	{
		if (!m_haveTrackingPreview)
		{
			//create the QImagePreviewReciever
			m_pTrackerPreview = new QTrackerPreviewReciever(m_parentWidget, m_targetLayout, QString::fromStdString(m_name));

			qRegisterMetaType<shared_ptr<TrackerDataSet>>("std::shared_ptr<TrackerDataSet>");
			connect(this, SIGNAL(previewReadyTracking(const std::shared_ptr<TrackerDataSet>)), m_pTrackerPreview, SLOT(previewReadyTracking(const std::shared_ptr<TrackerDataSet>)));
			m_haveTrackingPreview = true;
		}
	}

	std::tuple<double, double, bool> previewBuilderQT::computeWorldSize(std::shared_ptr <USImage> image)
	{
		double worldWidth;
		double worldHeight;
		bool keepRatio;
		if (image->getImageProperties()->getImageState() == USImageProperties::ImageState::Scan)
		{
			double resolution = image->getImageProperties()->getImageResolution();
			worldWidth = image->getSize().x * resolution;
			worldHeight = image->getSize().y * resolution;
			keepRatio = true;
		}
		else {
			worldWidth = image->getImageProperties()->getNumScanlines() - 1;
			worldHeight = image->getImageProperties()->getDepth();
			keepRatio = false;
		}
		return make_tuple(worldWidth, worldHeight, keepRatio);
	}

	tbb::flow::graph_node *
		previewBuilderQT::getInput(size_t index) {
		if (index == 0)
		{
			return &m_nodeIn;
		}
		else
		{
			return nullptr;
		}
	};
}