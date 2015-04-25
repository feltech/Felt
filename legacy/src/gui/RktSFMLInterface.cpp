#include <boost/format.hpp>
#include "gui/RktSFMLInterface.hpp"

using namespace felt;


RktSFMLInterface::~RktSFMLInterface()
{
	fprintf(stdout, "RktSFMLInterface::~RktSFMLInterface\n");
	try
	{
//		m_prktContext->RemoveReference();
		fprintf(stdout, "Shutdown\n");
		Rocket::Core::Shutdown();
	}
	catch (std::exception& e)
	{
		 fprintf(stderr, "Window::~Window error: %s\n", e.what());
		 fprintf(stdout, "Window::~Window error: %s\n", e.what());
		 throw e;
	}

	fprintf(stdout, "Window::~Window - done\n");
}


RktSFMLInterface::RktSFMLInterface(sf::RenderWindow* psfWindow) :
	m_psfWindow(NULL),
	m_prktContext(NULL),
	m_rktRenderer(),
	m_rktSystemInterface(),
	m_rktFileInterface("assets/")
{
	fprintf(stdout, "RktSFMLInterface::RktSFMLInterface v1\n");
	if (psfWindow == NULL)
	{
		fprintf(stdout, "psfWindow == NULL\n");
		throw std::runtime_error(
			str(boost::format(
				"RktSFMLInterface::RktSFMLInterface failed: %1%"
			) % psfWindow)
		);
	}

	m_psfWindow = psfWindow;
	m_rktRenderer.SetWindow(m_psfWindow);
	Rocket::Core::SetFileInterface(&m_rktFileInterface);
	Rocket::Core::SetRenderInterface(&m_rktRenderer);
	Rocket::Core::SetSystemInterface(&m_rktSystemInterface);
	if(!Rocket::Core::Initialise())
		throw std::runtime_error("Rocket::Core::Initialise failed");
}


void RktSFMLInterface::rktInit(std::string sfContextName)
{
	fprintf(stdout, "RktSFMLInterface::rktInit\n");
	m_prktContext = Rocket::Core::GetContext(sfContextName.c_str());

	if (m_prktContext == NULL)
	{
		fprintf(stdout, "m_prktContext == NULL\n");
		throw std::runtime_error(
			str(boost::format(
				"Rocket::Core::GetContext failed: failed to find \"%1%\""
			) % sfContextName)
		);
	}

	Rocket::Core::FontDatabase::LoadFontFace("Delicious-Bold.otf");
	Rocket::Core::FontDatabase::LoadFontFace("Delicious-BoldItalic.otf");
	Rocket::Core::FontDatabase::LoadFontFace("Delicious-Italic.otf");
	Rocket::Core::FontDatabase::LoadFontFace("Delicious-Roman.otf");

//	m_prktContext = Rocket::Core::CreateContext(
//		"main",
//		Rocket::Core::Vector2i(
//			m_psfWindow->getSize().x, m_psfWindow->getSize().y
//		)
//	);

}

void RktSFMLInterface::rktInitDebugger()
{
	Rocket::Debugger::Initialise(m_prktContext);
}

void RktSFMLInterface::rktEvent(sf::Event* sfEvent) {

	switch(sfEvent->type)
	{
	case sf::Event::Resized:
		//fprintf(stdout, "sf::Event::Resized\n");
		m_rktRenderer.Resize();
		break;
	case sf::Event::MouseMoved:
		//fprintf(stdout, "sf::Event::MouseMoved\n");
		m_prktContext->ProcessMouseMove(sfEvent->mouseMove.x, sfEvent->mouseMove.y,
			m_rktSystemInterface.GetKeyModifiers(m_psfWindow));
		break;
	case sf::Event::MouseButtonPressed:
		//fprintf(stdout, "sf::Event::MouseButtonPressed\n");
		m_prktContext->ProcessMouseButtonDown(sfEvent->mouseButton.button,
			m_rktSystemInterface.GetKeyModifiers(m_psfWindow));
		break;
	case sf::Event::MouseButtonReleased:
		//fprintf(stdout, "sf::Event::MouseButtonReleased\n");
		m_prktContext->ProcessMouseButtonUp(sfEvent->mouseButton.button,
			m_rktSystemInterface.GetKeyModifiers(m_psfWindow));
		break;
	case sf::Event::MouseWheelMoved:
		//fprintf(stdout, "sf::Event::MouseWheelMoved\n");
		m_prktContext->ProcessMouseWheel(-sfEvent->mouseWheel.delta,
			m_rktSystemInterface.GetKeyModifiers(m_psfWindow));
		break;
	case sf::Event::TextEntered:
		//fprintf(stdout, "sf::Event::TextEntered\n");
		if (sfEvent->text.unicode > 32)
			m_prktContext->ProcessTextInput(sfEvent->text.unicode);
		break;
	case sf::Event::KeyPressed:
		//fprintf(stdout, "sf::Event::KeyPressed\n");
		m_prktContext->ProcessKeyDown(m_rktSystemInterface.TranslateKey(sfEvent->key.code),
			m_rktSystemInterface.GetKeyModifiers(m_psfWindow));
		break;
	case sf::Event::KeyReleased:
		//fprintf(stdout, "sf::Event::KeyReleased\n");
		if(sfEvent->key.code == sf::Keyboard::F8)
		{
			Rocket::Debugger::SetVisible(!Rocket::Debugger::IsVisible());
		};

		if(sfEvent->key.code == sf::Keyboard::Escape) {
			m_psfWindow->close();
		}

		m_prktContext->ProcessKeyUp(m_rktSystemInterface.TranslateKey(sfEvent->key.code),
			m_rktSystemInterface.GetKeyModifiers(m_psfWindow));
		break;
	case sf::Event::Closed:
		//fprintf(stdout, "sf::Event::Closed\n");
		m_psfWindow->close();
		break;
	default:
		break;
	};
}
