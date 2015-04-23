/*
 * Window.hpp
 *
 *  Created on: 11 Apr 2015
 *      Author: dave
 */

#ifndef GUI_WINDOW_HPP_
#define GUI_WINDOW_HPP_

#include <GL/glew.h>
#include "Rocket/Core.h"
#include "gui/RktSFMLSystem.hpp"
#include "gui/RktSFMLRenderer.hpp"
#include "Rocket/Core/Input.h"
#include "Rocket/Debugger/Debugger.h"
#include "gui/RktFileInterface.hpp"

namespace felt
{
	class RktSFMLInterface
	{
	private:
		sf::RenderWindow* 			m_psfWindow;

		RktSFMLRenderer				m_rktRenderer;
		RktSFMLSystem				m_rktSystemInterface;
		RktFileInterface 			m_rktFileInterface;
		Rocket::Core::Context*		m_prktContext;

	public:
		RktSFMLInterface(sf::RenderWindow* psfWindow);
		~RktSFMLInterface();
		void 	rktInit(std::string sfContextName = "main");
		void 	rktEvent(sf::Event* sfEvent);
		void	rktInitDebugger();
	};
}
#endif /* GUI_WINDOW_HPP_ */
