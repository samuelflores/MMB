#ifndef _CONSTEXPR_H
#define _CONSTEXPR_H

#ifdef USE_MMB_CONSTEXPR
    #define MMB_CONSTEXPR constexpr
#else
    #define MMB_CONSTEXPR
#endif // USE_MMB_CONSTEXPR

#endif // _CONSTEXPR_H
