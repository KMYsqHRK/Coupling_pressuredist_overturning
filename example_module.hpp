#pragma once
#if !defined(PW_SDK_EXAMPLE_MODULE_H_INCLUDED)
#define PW_SDK_EXAMPLE_MODULE_H_INCLUDED

#include <memory>
#include "sdk/particleworks_api.hpp"

namespace pw { namespace example {

class Module {
protected:
	using UserFunctionPtr = std::shared_ptr<api::UserFunction>;

protected:
	api::Session						m_session;
	std::unordered_set<UserFunctionPtr>	m_functions;

public:
	virtual ~Module();

	virtual void initialize() = 0;

	virtual void run(const std::vector<std::string> & args);

protected:
	void add_user_function(const char * kernel, PW_CALL_POINT_t point, UserFunctionPtr uf);

public:
	static std::shared_ptr<Module> create(std::string key, PW_SESSION_t session);
	void set_session(PW_SESSION_t session){m_session=session;}
};

template <typename F>
inline PW_ERROR_code_t
safe_call(F && fn)
{
	try {
		fn();
	} catch (api::Error & e) {
		api::Logger(PW_LOG_warning_c) << "api error : " << e.what();
		return e.code();
	} catch (std::exception & e) {
		api::Logger(PW_LOG_warning_c) << "runtime error : " << e.what();
		return PW_ERROR_unknown;
	} catch (...) {
		api::Logger(PW_LOG_warning_c) << "unknown error";
		return PW_ERROR_unknown;
	}
	return PW_ERROR_no_error;
}

}}

#endif
