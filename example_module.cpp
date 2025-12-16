#include <memory>
#include "example_module.hpp"

namespace pw { namespace example {

Module::~Module()
{
	for (auto && uf : m_functions) {
		m_session.delete_user_function(uf.get());
	}
}

void
Module::add_user_function(const char * kernel, PW_CALL_POINT_t point, UserFunctionPtr uf)
{
	m_session.add_user_function(kernel, point, uf.get());
	m_functions.insert(uf);
	api::Logger(PW_LOG_info_c) << "User function installed : " << uf->name();
}

void
Module::run(const std::vector<std::string> & args)
{
	m_session.solver().step(-1);
}

}}
