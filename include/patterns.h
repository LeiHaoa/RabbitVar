#ifndef _PATTERNS_H
#define _PATTERNS_H

#include <regex>
#include <string>

class Patterns{
public:
	std::regex BEGIN_DIGITS{"^(\\d+)"};
    //Variation patterns
    //#define BEGIN_DIGITS  "^(\\d+)"
	std::regex MC_Z_NUM_S_ANY_NUM_S{"\\d+S\\S*\\d+S"};
	std::regex UP_NUMBER_END{"\\^(\\d+)$"};
	std::regex BEGIN_MINUS_NUMBER_ANY{"^-\\d+(.*)"};
	std::regex BEGIN_MINUS_NUMBER_CARET{"^-\\d+\\^"};
	//#define BEGIN_MINUS_NUMBER "^-(\\d+)"
	std::regex BEGIN_MINUS_NUMBER{"^-(\\d+).*"};
    std::regex MINUS_NUM_NUM{"-\\d\\d"};
    std::regex HASH_GROUP_CARET_GROUP{"#(.+)\\^(.+)"};

    //Sclip patterns
    std::regex B_A7{"^.AAAAAAA.*"};
    std::regex B_T7{"^.TTTTTTT.*"};

    //ATGC patterns
    std::regex CARET_ATGNC{"\\^([ATGNC]+)"};
    std::regex CARET_ATGC_END {"\\^([ATGC]+)$"};
    std::regex AMP_ATGC {".*&([ATGC]+).*"};
    std::regex BEGIN_PLUS_ATGC{"^\\+([ATGC]+).*"};
    std::regex HASH_ATGC{"#([ATGC]+).*"};
    std::regex ATGSs_AMP_ATGSs_END{"(\\+[ATGC]+)&[ATGC]+$"};
    std::regex MINUS_NUMBER_AMP_ATGCs_END {"(-\\d+)&[ATGC]+$"};
    std::regex MINUS_NUMBER_ATGNC_SV_ATGNC_END{"^-\\d+\\^([ATGNC]+)<...\\d+>([ATGNC]+)$"};
    std::regex BEGIN_ATGC_END{"^[ATGC]+$"};

    //SV patterns
    std::regex DUP_NUM{"<dup(\\d+)"};
    std::regex DUP_NUM_ATGC{"<dup(\\d+)>([ATGC]+)$"};
    std::regex INV_NUM{"<inv(\\d+)"};
    std::regex SOME_SV_NUMBERS{"<(...)\\d+>"};
    std::regex ANY_SV{"<(...)>"};

    //File and columns patterns
    std::regex SAMPLE_PATTERN{"([^\\/\\._]+).sorted[^\\/]*.bam"};
    std::regex SAMPLE_PATTERN2{"([^\\/]+)[_\\.][^\\/]*bam"};
    std::regex INTEGER_ONLY{"^\\d+$"};

    //CIGAR patterns
    /**
     * Std::Regexp finds number followed by S followed by number and I or D at the start of string
     */
    std::regex BEGIN_NUMBER_S_NUMBER_IorD{"^(\\d+)S(\\d+)([ID])"};
    /**
     * Std::Regexp finds number followed by I or D followed by number followed by S followed by end of string
     */
    std::regex NUMBER_IorD_NUMBER_S_END{"(\\d+)([ID])(\\d+)S$"};
    /**
     * Std::Regexp finds number
     * followed by S followed by number followed by M followed by number followed by I or D at the start of string
     */
    std::regex BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD{"^(\\d+)S(\\d+)M(\\d+)([ID])"};
    /**
     * Std::Regexp finds number
     * followed by I or D followed by number followed by M followed by number followed by S followed by end of string
     */
    std::regex NUMBER_IorD_NUMBER_M_NUMBER_S_END{"(\\d+)([ID])(\\d+)M(\\d+)S$"};
    /**
     * Std::Regexp finds digit-M-number-(I or D)-number-M
     */
    std::regex BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M{"^(\\d)M(\\d+)([ID])(\\d+)M"};

    std::regex BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_{"^\\dM\\d+[ID]\\d+M"};
    /**
     * Std::Regexp replaces number-(I or D)-number-M
     */
    std::regex NUMBER_IorD_DIGIT_M_END{"(\\d+)([ID])(\\d)M$"};
    std::regex NUMBER_IorD_NUMBER_M_END{"(\\d+)([ID])(\\d+)M$"};
    /**
     * Std::Regexp finds number-M-number-D-digit-M-number-I-digit-M-number-D. ^(.*?) captures everything before the
     * M-D-M-I-M-D complex
     */
    std::regex D_M_D_DD_M_D_I_D_M_D_DD{"^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M"};
    std::regex D_M_D_DD_M_D_I_D_M_D_DD_prim{"(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M"};
    //------------[haoz:] ??? -------------
    //std::regex threeDeletionsPattern{"^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)" + "M"};
    std::regex threeDeletionsPattern{"^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M"};
    std::regex threeIndelsPattern{"^(.*?)(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M"};
    std::regex DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM{"\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M"};
    std::regex DM_DD_DM_DD_DM_DD_DM{"\\d+M\\d+D\\d+M\\d+D\\d+M\\d+D\\d+M"};
    /**
     * Std::Regexp finds number-D-digit-M-number-D and optionally number-I
     */
    std::regex DIG_D_DIG_M_DIG_DI_DIGI{"(\\d+)D(\\d+)M(\\d+)([DI])(\\d+I)?"};
    /**
     * Std::Regexp finds number-I-digit-M-number and optionally number-I
     */
    std::regex DIG_I_DIG_M_DIG_DI_DIGI{"(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?"};
    /**
     * Std::Regexp finds not_digit-number-I-number-M-number-D and optionally number-I
     */
    std::regex NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI{"(\\D)(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?"};
    /**
     * Std::Regexp finds number-D-number-D
     */
    std::regex DIG_D_DIG_D{"(\\d+)D(\\d+)D"};
    /**
     * Std::Regexp finds number-I-number-I
     */
    std::regex DIG_I_DIG_I{"(\\d+)I(\\d+)I"};
    std::regex BEGIN_ANY_DIG_M_END{"^(.*?)(\\d+)M$"};

    std::regex DIG_M_END{"\\d+M$"};
    std::regex BEGIN_DIG_M{"^(\\d+)M"};
    /**
     * Std::Regexp finds number-S-number-M at start of CIGAR string
     */
    std::regex DIG_S_DIG_M{"^(\\d+)S(\\d+)M"};
    std::regex DIG_M_DIG_S_END{"\\d+M\\d+S$"};
    /**
     * Std::Regexp finds number-M-number-S at end of CIGAR string ^(.*?) captures everything before M-S complex
     */
    std::regex ANY_NUMBER_M_NUMBER_S_END{"^(.*?)(\\d+)M(\\d+)S$"};
    /**
     * Std::Regexp finds number followed by D at the start of string
     */
    std::regex BEGIN_NUMBER_D{"^(\\d+)D"};
    std::regex END_NUMBER_D{"(\\d+)D$"};
    /**
     * Std::Regexp finds number followed by I at the start of string
     */
    std::regex BEGIN_NUMBER_I{"^(\\d+)I"};
    /**
     * Std::Regexp finds number followed by I at the end of string
     */
    std::regex END_NUMBER_I{"(\\d+)I$"};

    /**
     * Std::Regexp finds numbers followed by M (matched) or D (deleted) in CIGAR string
     */
    std::regex ALIGNED_LENGTH_MND{"(\\d+)[MND]"};

    /**
     * The total aligned length, excluding soft-clipped bases and insertions
     */
    std::regex ALIGNED_LENGTH_MD{"(\\d+)[MD=X]"};

    std::regex SOFT_CLIPPED{"(\\d+)[MIS]"};
    std::regex SA_CIGAR_D_S_5clip{"^\\d\\d+S"};
    std::regex SA_CIGAR_D_S_3clip{"\\d\\dS$"};
    std::regex BEGIN_dig_dig_S_ANY_dig_dig_S_END{"^\\d\\dS.*\\d\\dS$"};
    std::regex BEGIN_NUM_S_OR_BEGIN_NUM_H{"^(\\d+)S|^\\d+H"};
	std::regex END_NUM_S_OR_NUM_H{"(\\d+)S$|H$"};
};
#endif
