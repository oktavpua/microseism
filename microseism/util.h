#ifndef UTIL_H
#define UTIL_H

/** 
 * \file   util.h
 * \author Денис Демидов <ddemidov@ksu.ru>
 * \brief  Вспомогательные функции.
 */

#include <iostream>
#include <string>
#include <stdexcept>

/// Логгер.
class Log {
    public:
	/// Признак окончания строки.
	static const struct endl_tag{} endl;

	/// Признак очистки буфера.
	static const struct flush_tag{} flush;

	/// Конструктор.
	/**
	 * \param os   Поток, в который будет выводиться лог.
	 * \param rank Номер процесса.
	 *
	 * \note Лог выводится только для нулевого процесса.
	 */
	Log(std::ostream &os, int rank);

	/// Вывод значения в лог.
	/**
	 * \param msg Значение.
	 */
	template <class T>
	Log& operator<<(T &&msg) {
	    if (rank == 0)
		os << std::forward<T>(msg);
	    return *this;
	}

	Log& operator<<(std::ios_base& manip(std::ios_base&)) {
	    if (rank == 0)
		os << manip;
	    return *this;
	}

	/// Перенос строки.
	Log& operator<<(endl_tag  endl);

	/// Очистка буфера.
	Log& operator<<(flush_tag flush);
    private:
	std::ostream &os;
	int           rank;
};

/// Проверка условия на истинность.
/**
 * В вслучае, если условие не выполняется, выбрасывается исключение с
 * соответствующим сообщением об ошибке.
 *
 * \param cond     Условие.
 * \param fail_msg Сообщение об ошибке.
 */
template <class Cond>
inline void precondition(const Cond &cond, const std::string &fail_msg) {
    if (!static_cast<bool>(cond)) throw std::logic_error(fail_msg);
}

#endif
