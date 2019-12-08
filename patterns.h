#ifndef _PATTERNS_H
#define _PATTERNS_H

#include <regex>
#include <string>

class Patterns{
public:
	regex BEGIN_DIGITS{"^(\\d+)"};
    //Variation patterns
    //#define BEGIN_DIGITS  "^(\\d+)"
	regex MC_Z_NUM_S_ANY_NUM_S{"\\d+S\\S*\\d+S"};
	regex UP_NUMBER_END{"\\^(\\d+)$"};
	regex BEGIN_MINUS_NUMBER_ANY{"^-\\d+(.*)"};
	regex BEGIN_MINUS_NUMBER_CARET{"^-\\d+\\^"};
	//#define BEGIN_MINUS_NUMBER "^-(\\d+)"
	regex BEGIN_MINUS_NUMBER{"^-(\\d+).*"};
    regex MINUS_NUM_NUM{"-\\d\\d"};
    regex HASH_GROUP_CARET_GROUP{"#(.+)\\^(.+)"};

    //Sclip patterns
    regex B_A7{"^.AAAAAAA.*"};
    regex B_T7{"^.TTTTTTT.*"};

    //ATGC patterns
    regex CARET_ATGNC{"\\^([ATGNC]+)"};
    regex CARET_ATGC_END {"\\^([ATGC]+)$"};
    regex AMP_ATGC {".*&([ATGC]+).*"};
    regex BEGIN_PLUS_ATGC{"^\\+([ATGC]+).*"};
    regex HASH_ATGC{"#([ATGC]+).*"};
    regex ATGSs_AMP_ATGSs_END{"(\\+[ATGC]+)&[ATGC]+$"};
    regex MINUS_NUMBER_AMP_ATGCs_END {"(-\\d+)&[ATGC]+$"};
    regex MINUS_NUMBER_ATGNC_SV_ATGNC_END{"^-\\d+\\^([ATGNC]+)<...\\d+>([ATGNC]+)$"};
    regex BEGIN_ATGC_END{"^[ATGC]+$"};

    //SV patterns
    regex DUP_NUM{"<dup(\\d+)"};
    regex DUP_NUM_ATGC{"<dup(\\d+)>([ATGC]+)$"};
    regex INV_NUM{"<inv(\\d+)"};
    regex SOME_SV_NUMBERS{"<(...)\\d+>"};
    regex ANY_SV{"<(...)>"};

    //File and columns patterns
    regex SAMPLE_PATTERN{"([^\\/\\._]+).sorted[^\\/]*.bam"};
    regex SAMPLE_PATTERN2{"([^\\/]+)[_\\.][^\\/]*bam"};
    regex INTEGER_ONLY{"^\\d+$"};

    //CIGAR patterns
    /**
     * Regexp finds number followed by S followed by number and I or D at the start of string
     */
    regex BEGIN_NUMBER_S_NUMBER_IorD{"^(\\d+)S(\\d+)([ID])"};
    /**
     * Regexp finds number followed by I or D followed by number followed by S followed by end of string
     */
    regex NUMBER_IorD_NUMBER_S_END{"(\\d+)([ID])(\\d+)S$"};
    /**
     * Regexp finds number
     * followed by S followed by number followed by M followed by number followed by I or D at the start of string
     */
    regex BEGIN_NUMBER_S_NUMBER_M_NUMBER_IorD{"^(\\d+)S(\\d+)M(\\d+)([ID])"};
    /**
     * Regexp finds number
     * followed by I or D followed by number followed by M followed by number followed by S followed by end of string
     */
    regex NUMBER_IorD_NUMBER_M_NUMBER_S_END{"(\\d+)([ID])(\\d+)M(\\d+)S$"};
    /**
     * Regexp finds digit-M-number-(I or D)-number-M
     */
    regex BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M{"^(\\d)M(\\d+)([ID])(\\d+)M"};

    regex BEGIN_DIGIT_M_NUMBER_IorD_NUMBER_M_{"^\\dM\\d+[ID]\\d+M"};
    /**
     * Regexp replaces number-(I or D)-number-M
     */
    regex NUMBER_IorD_DIGIT_M_END{"(\\d+)([ID])(\\d)M$"};
    regex NUMBER_IorD_NUMBER_M_END{"(\\d+)([ID])(\\d+)M$"};
    /**
     * Regexp finds number-M-number-D-digit-M-number-I-digit-M-number-D. ^(.*?) captures everything before the
     * M-D-M-I-M-D complex
     */
    regex D_M_D_DD_M_D_I_D_M_D_DD{"^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M"};
    regex D_M_D_DD_M_D_I_D_M_D_DD_prim{"(\\d+)M(\\d+)D(\\d+)M(\\d+)I(\\d+)M(\\d+)D(\\d+)M"};
    //------------[haoz:] ??? -------------
    //regex threeDeletionsPattern{"^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)" + "M"};
    regex threeDeletionsPattern{"^(.*?)(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M(\\d+)D(\\d+)M"};
    regex threeIndelsPattern{"^(.*?)(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M(\\d+)([DI])(\\d+)M"};
    regex DIGM_D_DI_DIGM_D_DI_DIGM_DI_DIGM{"\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M\\d+[DI]\\d+M"};
    regex DM_DD_DM_DD_DM_DD_DM{"\\d+M\\d+D\\d+M\\d+D\\d+M\\d+D\\d+M"};
    /**
     * Regexp finds number-D-digit-M-number-D and optionally number-I
     */
    regex DIG_D_DIG_M_DIG_DI_DIGI{"(\\d+)D(\\d+)M(\\d+)([DI])(\\d+I)?"};
    /**
     * Regexp finds number-I-digit-M-number and optionally number-I
     */
    regex DIG_I_DIG_M_DIG_DI_DIGI{"(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?"};
    /**
     * Regexp finds not_digit-number-I-number-M-number-D and optionally number-I
     */
    regex NOTDIG_DIG_I_DIG_M_DIG_DI_DIGI{"(\\D)(\\d+)I(\\d+)M(\\d+)([DI])(\\d+I)?"};
    /**
     * Regexp finds number-D-number-D
     */
    regex DIG_D_DIG_D{"(\\d+)D(\\d+)D"};
    /**
     * Regexp finds number-I-number-I
     */
    regex DIG_I_DIG_I{"(\\d+)I(\\d+)I"};
    regex BEGIN_ANY_DIG_M_END{"^(.*?)(\\d+)M$"};

    regex DIG_M_END{"\\d+M$"};
    regex BEGIN_DIG_M{"^(\\d+)M"};
    /**
     * Regexp finds number-S-number-M at start of CIGAR string
     */
    regex DIG_S_DIG_M{"^(\\d+)S(\\d+)M"};
    regex DIG_M_DIG_S_END{"\\d+M\\d+S$"};
    /**
     * Regexp finds number-M-number-S at end of CIGAR string ^(.*?) captures everything before M-S complex
     */
    regex ANY_NUMBER_M_NUMBER_S_END{"^(.*?)(\\d+)M(\\d+)S$"};
    /**
     * Regexp finds number followed by D at the start of string
     */
    regex BEGIN_NUMBER_D{"^(\\d+)D"};
    regex END_NUMBER_D{"(\\d+)D$"};
    /**
     * Regexp finds number followed by I at the start of string
     */
    regex BEGIN_NUMBER_I{"^(\\d+)I"};
    /**
     * Regexp finds number followed by I at the end of string
     */
    regex END_NUMBER_I{"(\\d+)I$"};

    /**
     * Regexp finds numbers followed by M (matched) or D (deleted) in CIGAR string
     */
    regex ALIGNED_LENGTH_MND{"(\\d+)[MND]"};

    /**
     * The total aligned length, excluding soft-clipped bases and insertions
     */
    regex ALIGNED_LENGTH_MD{"(\\d+)[MD=X]"};

    regex SOFT_CLIPPED{"(\\d+)[MIS]"};
    regex SA_CIGAR_D_S_5clip{"^\\d\\d+S"};
    regex SA_CIGAR_D_S_3clip{"\\d\\dS$"};
    regex BEGIN_dig_dig_S_ANY_dig_dig_S_END{"^\\d\\dS.*\\d\\dS$"};
    regex BEGIN_NUM_S_OR_BEGIN_NUM_H{"^(\\d+)S|^\\d+H"};
	regex END_NUM_S_OR_NUM_H{"(\\d+)S$|H$"};
};
#endif
