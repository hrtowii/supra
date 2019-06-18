// ================================================================================================
// 
// Copyright (C) 2017, Rüdiger Göbl - all rights reserved
// Copyright (c) 2019, NVIDIA CORPORATION. All rights reserved.
//
//          Rüdiger Göbl
//          Email r.goebl@tum.de
//          Chair for Computer Aided Medical Procedures
//          Technische Universität München
//          Boltzmannstr. 3, 85748 Garching b. München, Germany
// 
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License, version 2.1, as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this program.  If not, see
// <http://www.gnu.org/licenses/>.
//
// ================================================================================================

#ifndef __BEAMFORMINGMVPCGNODE_H__
#define __BEAMFORMINGMVPCGNODE_H__

#ifdef HAVE_BEAMFORMER_MINIMUM_VARIANCE

#include <memory>
#include <tbb/flow_graph.h>

#include "AbstractNode.h"
#include "RecordObject.h"

namespace supra
{
	//forward declarations
	class USImageProperties;
	class USImage;
	class USRawData;

	class BeamformingMVpcgNode : public AbstractNode {
	public:
		BeamformingMVpcgNode(tbb::flow::graph& graph, const std::string & nodeID, bool queueing);
		~BeamformingMVpcgNode();

		virtual size_t getNumInputs() { return 1; }
		virtual size_t getNumOutputs() { return 1; }

		virtual tbb::flow::graph_node * getInput(size_t index) {
			if (index == 0)
			{
				return m_node.get();
			}
			return nullptr;
		};

		virtual tbb::flow::graph_node * getOutput(size_t index) {
			if (index == 0)
			{
				return m_node.get();
			}
			return nullptr;
		};

	protected:
		void configurationChanged();
		void configurationEntryChanged(const std::string& configKey);

	private:
		std::shared_ptr<RecordObject> checkTypeAndBeamform(std::shared_ptr<RecordObject> mainObj);
		template <typename RawDataType>
		std::shared_ptr<USImage> beamformTemplated(std::shared_ptr<const USRawData> rawData);
		void updateImageProperties(std::shared_ptr<const USImageProperties> imageProperties);

		std::shared_ptr<const USImageProperties> m_lastSeenImageProperties;
		std::shared_ptr<USImageProperties> m_editedImageProperties;

		std::mutex m_mutex;

		std::unique_ptr<tbb::flow::graph_node> m_node;

		uint32_t m_subArraySize;
		uint32_t m_temporalSmoothing;
		uint32_t m_maxIterations;
		double m_convergenceThreshold;
		double m_outputClamp;
		DataType m_outputType;
		uint32_t m_maxIterationsOverride;
		double m_subArrayScalingPower;
		bool m_computeMeans;
	};
}

#endif //HAVE_BEAMFORMER_MINIMUM_VARIANCE

#endif //!__BEAMFORMINGMVPCGNODE_H__
