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

#ifndef __TRACKERINTERFACESIMULATED_H__
#define __TRACKERINTERFACESIMULATED_H__

#ifdef HAVE_DEVICE_TRACKING_SIM

#include <AbstractInput.h>
#include <TrackerDataSet.h>

#include <mutex>

namespace supra
{
	class TrackerInterfaceSimulated : public AbstractInput
	{
	public:

		TrackerInterfaceSimulated(tbb::flow::graph & graph, const std::string & nodeID)
			: AbstractInput(graph, nodeID,1)
			, m_frozen(false)
		{
			m_valueRangeDictionary.set<int>("frequency", 5, 200, 50, "Frequency");
			m_valueRangeDictionary.set<int>("trackerID", 123, "Tacker ID");
		};

		~TrackerInterfaceSimulated() {};

		//Functions to be overwritten
	public:
		virtual void initializeDevice() { };
		virtual bool ready() { return true; };

		virtual std::vector<size_t> getImageOutputPorts() { return{}; };
		virtual std::vector<size_t> getTrackingOutputPorts() { return{ 0 }; };

		virtual void freeze();
		virtual void unfreeze();
	protected:
		virtual void startAcquisition();
		//Needs to be thread safe
		virtual void stopAcquisition() {};
		//Needs to be thread safe
		virtual void configurationEntryChanged(const std::string& configKey);
		//Needs to be thread safe
		virtual void configurationChanged();

		virtual bool timerCallback();

	private:
		double m_initTime; 
		int m_frequency;
		std::mutex m_objectMutex;
		std::atomic_bool m_frozen;
		int m_trackerId;
		

		void readConfiguration();
	};
}

#endif //!HAVE_DEVICE_TRACKING_SIM

#endif //!__TRACKERINTERFACESIMULATED_H__
