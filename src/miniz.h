#include "quakedef.h"

#define Q_MINIZ_H 1
#define MINIZ_EXPORT
#define MINIZ_NO_STDIO
#define MINIZ_NO_TIME
#define MINIZ_NO_DEFLATE_APIS
#define MINIZ_NO_ARCHIVE_WRITING_APIS
#define MINIZ_NO_ZLIB_APIS
#define MINIZ_NO_ZLIB_COMPATIBLE_NAMES
#define MINIZ_LITTLE_ENDIAN (SDL_BYTEORDER == SDL_LIL_ENDIAN)
#define MINIZ_NO_ARCHIVE_WRITING_APIS
#define MINIZ_USE_UNALIGNED_LOADS_AND_STORES 1
#define MINIZ_UNALIGNED_USE_MEMCPY
#define MINIZ_HAS_64BIT_REGISTERS 1
#define MZ_CRC32_INIT (0)
#define MZ_DEFLATED 8
#define MZ_VERSION "10.2.0"
#define MZ_VERNUM 0xA100
#define MZ_VER_MAJOR 10
#define MZ_VER_MINOR 2
#define MZ_VER_REVISION 0
#define MZ_VER_SUBREVISION 0
#define MZ_FALSE (0)
#define MZ_TRUE (1)
#define MZ_MACRO_END while (0)
#define MZ_FILE void *
#define MZ_ASSERT(x) assert(x)
#define MZ_MALLOC(x) malloc(x)
#define MZ_FREE(x) free(x)
#define MZ_REALLOC(p, x) realloc(p, x)
#define MZ_MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MZ_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MZ_CLEAR_OBJ(obj) memset(&(obj), 0, sizeof(obj))
#define MZ_CLEAR_ARR(obj) memset((obj), 0, sizeof(obj))
#define MZ_CLEAR_PTR(obj) memset((obj), 0, sizeof(*obj))
#define MZ_READ_LE16(p) *((const mz_uint16 *)(p))
#define MZ_READ_LE32(p) *((const mz_uint32 *)(p))
#define MZ_READ_LE64(p) (((mz_uint64)MZ_READ_LE32(p)) | (((mz_uint64)MZ_READ_LE32((const mz_uint8 *)(p) + sizeof(mz_uint32))) << 32U))
#define MZ_FORCEINLINE inline
#define MZ_UINT16_MAX (0xFFFFU)
#define MZ_UINT32_MAX (0xFFFFFFFFU)
#define TINFL_LZ_DICT_SIZE 32768
// Initializes the decompressor to its initial state.
#define tinfl_init(r)     \
    do                    \
    {                     \
        (r)->m_state = 0; \
    }                     \
    MZ_MACRO_END
#define tinfl_get_adler32(r) (r)->m_check_adler32
#define TINFL_USE_64BIT_BITBUF 1
#define TINFL_BITBUF_SIZE (64)

typedef unsigned long mz_ulong;
typedef void *(*mz_alloc_func)(void *opaque, size_t items, size_t size);
typedef void (*mz_free_func)(void *opaque, void *address);
typedef void *(*mz_realloc_func)(void *opaque, void *address, size_t items, size_t size);
typedef unsigned char mz_uint8;
typedef signed short mz_int16;
typedef unsigned short mz_uint16;
typedef unsigned int mz_uint32;
typedef unsigned int mz_uint;
typedef int64_t mz_int64;
typedef uint64_t mz_uint64;
typedef int mz_bool;
enum {
    TINFL_FLAG_PARSE_ZLIB_HEADER = 1,
    TINFL_FLAG_HAS_MORE_INPUT = 2,
    TINFL_FLAG_USING_NON_WRAPPING_OUTPUT_BUF = 4,
    TINFL_FLAG_COMPUTE_ADLER32 = 8
};
struct tinfl_decompressor_tag;
typedef struct tinfl_decompressor_tag tinfl_decompressor;
typedef enum {
    TINFL_STATUS_FAILED_CANNOT_MAKE_PROGRESS = -4,
    TINFL_STATUS_BAD_PARAM = -3,
    TINFL_STATUS_ADLER32_MISMATCH = -2,
    TINFL_STATUS_FAILED = -1,
    TINFL_STATUS_DONE = 0,
    TINFL_STATUS_NEEDS_MORE_INPUT = 1,
    TINFL_STATUS_HAS_MORE_OUTPUT = 2
} tinfl_status;
MINIZ_EXPORT tinfl_status tinfl_decompress(tinfl_decompressor *r, const mz_uint8 *pIn_buf_next, size_t *pIn_buf_size, mz_uint8 *pOut_buf_start, mz_uint8 *pOut_buf_next, size_t *pOut_buf_size, const mz_uint32 decomp_flags);
enum {
    TINFL_MAX_HUFF_TABLES = 3,
    TINFL_MAX_HUFF_SYMBOLS_0 = 288,
    TINFL_MAX_HUFF_SYMBOLS_1 = 32,
    TINFL_MAX_HUFF_SYMBOLS_2 = 19,
    TINFL_FAST_LOOKUP_BITS = 10,
    TINFL_FAST_LOOKUP_SIZE = 1 << TINFL_FAST_LOOKUP_BITS
};
typedef mz_uint64 tinfl_bit_buf_t;
struct tinfl_decompressor_tag {
    mz_uint32 m_state, m_num_bits, m_zhdr0, m_zhdr1, m_z_adler32, m_final, m_type, m_check_adler32, m_dist, m_counter, m_num_extra, m_table_sizes[TINFL_MAX_HUFF_TABLES];
    tinfl_bit_buf_t m_bit_buf;
    size_t m_dist_from_out_buf_start;
    mz_int16 m_look_up[TINFL_MAX_HUFF_TABLES][TINFL_FAST_LOOKUP_SIZE];
    mz_int16 m_tree_0[TINFL_MAX_HUFF_SYMBOLS_0 * 2];
    mz_int16 m_tree_1[TINFL_MAX_HUFF_SYMBOLS_1 * 2];
    mz_int16 m_tree_2[TINFL_MAX_HUFF_SYMBOLS_2 * 2];
    mz_uint8 m_code_size_0[TINFL_MAX_HUFF_SYMBOLS_0];
    mz_uint8 m_code_size_1[TINFL_MAX_HUFF_SYMBOLS_1];
    mz_uint8 m_code_size_2[TINFL_MAX_HUFF_SYMBOLS_2];
    mz_uint8 m_raw_header[4], m_len_codes[TINFL_MAX_HUFF_SYMBOLS_0 + TINFL_MAX_HUFF_SYMBOLS_1 + 137];
};
enum {
    MZ_ZIP_MAX_IO_BUF_SIZE = 64 * 1024,
    MZ_ZIP_MAX_ARCHIVE_FILENAME_SIZE = 512,
    MZ_ZIP_MAX_ARCHIVE_FILE_COMMENT_SIZE = 512
};
typedef struct {
    mz_uint32 m_file_index;
    mz_uint64 m_central_dir_ofs;
    mz_uint16 m_version_made_by;
    mz_uint16 m_version_needed;
    mz_uint16 m_bit_flag;
    mz_uint16 m_method;
    mz_uint32 m_crc32;
    mz_uint64 m_comp_size;
    mz_uint64 m_uncomp_size;
    mz_uint16 m_internal_attr;
    mz_uint32 m_external_attr;
    mz_uint64 m_local_header_ofs;
    mz_uint32 m_comment_size;
    mz_bool m_is_directory;
    mz_bool m_is_encrypted;
    mz_bool m_is_supported;
    char m_filename[MZ_ZIP_MAX_ARCHIVE_FILENAME_SIZE];
    char m_comment[MZ_ZIP_MAX_ARCHIVE_FILE_COMMENT_SIZE];
} mz_zip_archive_file_stat;
typedef size_t (*mz_file_read_func)(void *pOpaque, mz_uint64 file_ofs, void *pBuf, size_t n);
typedef size_t (*mz_file_write_func)(void *pOpaque, mz_uint64 file_ofs, const void *pBuf, size_t n);
typedef mz_bool (*mz_file_needs_keepalive)(void *pOpaque);
struct mz_zip_internal_state_tag;
typedef struct mz_zip_internal_state_tag mz_zip_internal_state;
typedef enum {
    MZ_ZIP_MODE_INVALID = 0,
    MZ_ZIP_MODE_READING = 1,
    MZ_ZIP_MODE_WRITING = 2,
    MZ_ZIP_MODE_WRITING_HAS_BEEN_FINALIZED = 3
} mz_zip_mode;
typedef enum {
    MZ_ZIP_FLAG_CASE_SENSITIVE = 0x0100,
    MZ_ZIP_FLAG_IGNORE_PATH = 0x0200,
    MZ_ZIP_FLAG_COMPRESSED_DATA = 0x0400,
    MZ_ZIP_FLAG_DO_NOT_SORT_CENTRAL_DIRECTORY = 0x0800,
    MZ_ZIP_FLAG_VALIDATE_LOCATE_FILE_FLAG = 0x1000,
    MZ_ZIP_FLAG_VALIDATE_HEADERS_ONLY = 0x2000,
    MZ_ZIP_FLAG_WRITE_ZIP64 = 0x4000,
    MZ_ZIP_FLAG_WRITE_ALLOW_READING = 0x8000,
    MZ_ZIP_FLAG_ASCII_FILENAME = 0x10000,
    MZ_ZIP_FLAG_WRITE_HEADER_SET_SIZE = 0x20000
} mz_zip_flags;
typedef enum {
    MZ_ZIP_TYPE_INVALID = 0,
    MZ_ZIP_TYPE_USER,
    MZ_ZIP_TYPE_MEMORY,
    MZ_ZIP_TYPE_HEAP,
    MZ_ZIP_TYPE_FILE,
    MZ_ZIP_TYPE_CFILE,
    MZ_ZIP_TOTAL_TYPES
} mz_zip_type;
typedef enum {
    MZ_ZIP_NO_ERROR = 0,
    MZ_ZIP_UNDEFINED_ERROR,
    MZ_ZIP_TOO_MANY_FILES,
    MZ_ZIP_FILE_TOO_LARGE,
    MZ_ZIP_UNSUPPORTED_METHOD,
    MZ_ZIP_UNSUPPORTED_ENCRYPTION,
    MZ_ZIP_UNSUPPORTED_FEATURE,
    MZ_ZIP_FAILED_FINDING_CENTRAL_DIR,
    MZ_ZIP_NOT_AN_ARCHIVE,
    MZ_ZIP_INVALID_HEADER_OR_CORRUPTED,
    MZ_ZIP_UNSUPPORTED_MULTIDISK,
    MZ_ZIP_DECOMPRESSION_FAILED,
    MZ_ZIP_COMPRESSION_FAILED,
    MZ_ZIP_UNEXPECTED_DECOMPRESSED_SIZE,
    MZ_ZIP_CRC_CHECK_FAILED,
    MZ_ZIP_UNSUPPORTED_CDIR_SIZE,
    MZ_ZIP_ALLOC_FAILED,
    MZ_ZIP_FILE_OPEN_FAILED,
    MZ_ZIP_FILE_CREATE_FAILED,
    MZ_ZIP_FILE_WRITE_FAILED,
    MZ_ZIP_FILE_READ_FAILED,
    MZ_ZIP_FILE_CLOSE_FAILED,
    MZ_ZIP_FILE_SEEK_FAILED,
    MZ_ZIP_FILE_STAT_FAILED,
    MZ_ZIP_INVALID_PARAMETER,
    MZ_ZIP_INVALID_FILENAME,
    MZ_ZIP_BUF_TOO_SMALL,
    MZ_ZIP_INTERNAL_ERROR,
    MZ_ZIP_FILE_NOT_FOUND,
    MZ_ZIP_ARCHIVE_TOO_LARGE,
    MZ_ZIP_VALIDATION_FAILED,
    MZ_ZIP_WRITE_CALLBACK_FAILED,
    MZ_ZIP_TOTAL_ERRORS
} mz_zip_error;
typedef struct {
    mz_uint64 m_archive_size;
    mz_uint64 m_central_directory_file_ofs;
    mz_uint32 m_total_files;
    mz_zip_mode m_zip_mode;
    mz_zip_type m_zip_type;
    mz_zip_error m_last_error;
    mz_uint64 m_file_offset_alignment;
    mz_alloc_func m_pAlloc;
    mz_free_func m_pFree;
    mz_realloc_func m_pRealloc;
    void *m_pAlloc_opaque;
    mz_file_read_func m_pRead;
    mz_file_write_func m_pWrite;
    mz_file_needs_keepalive m_pNeeds_keepalive;
    void *m_pIO_opaque;
    mz_zip_internal_state *m_pState;
} mz_zip_archive;
MINIZ_EXPORT mz_bool mz_zip_reader_init(mz_zip_archive *pZip, mz_uint64 size, mz_uint flags);
MINIZ_EXPORT mz_bool mz_zip_reader_end(mz_zip_archive *pZip);
MINIZ_EXPORT mz_zip_mode mz_zip_get_mode(mz_zip_archive *pZip);
MINIZ_EXPORT mz_zip_type mz_zip_get_type(mz_zip_archive *pZip);
MINIZ_EXPORT mz_uint mz_zip_reader_get_num_files(mz_zip_archive *pZip);
MINIZ_EXPORT mz_uint64 mz_zip_get_archive_size(mz_zip_archive *pZip);
MINIZ_EXPORT MZ_FILE *mz_zip_get_cfile(mz_zip_archive *pZip);
MINIZ_EXPORT size_t mz_zip_read_archive_data(mz_zip_archive *pZip, mz_uint64 file_ofs, void *pBuf, size_t n);
MINIZ_EXPORT mz_zip_error mz_zip_set_last_error(mz_zip_archive *pZip, mz_zip_error err_num);
MINIZ_EXPORT mz_zip_error mz_zip_peek_last_error(mz_zip_archive *pZip);
MINIZ_EXPORT mz_zip_error mz_zip_clear_last_error(mz_zip_archive *pZip);
MINIZ_EXPORT mz_zip_error mz_zip_get_last_error(mz_zip_archive *pZip);
MINIZ_EXPORT const char *mz_zip_get_error_string(mz_zip_error mz_err);
MINIZ_EXPORT mz_bool mz_zip_reader_is_file_a_directory(mz_zip_archive *pZip, mz_uint file_index);
MINIZ_EXPORT mz_bool mz_zip_reader_is_file_encrypted(mz_zip_archive *pZip, mz_uint file_index);
MINIZ_EXPORT mz_bool mz_zip_reader_is_file_supported(mz_zip_archive *pZip, mz_uint file_index);
MINIZ_EXPORT mz_uint mz_zip_reader_get_filename(mz_zip_archive *pZip, mz_uint file_index, char *pFilename, mz_uint filename_buf_size);
MINIZ_EXPORT int mz_zip_reader_locate_file(mz_zip_archive *pZip, const char *pName, const char *pComment, mz_uint flags);
MINIZ_EXPORT mz_bool mz_zip_reader_locate_file_v2(mz_zip_archive *pZip, const char *pName, const char *pComment, mz_uint flags, mz_uint32 *file_index);
MINIZ_EXPORT mz_bool mz_zip_reader_file_stat(mz_zip_archive *pZip, mz_uint file_index, mz_zip_archive_file_stat *pStat);
MINIZ_EXPORT mz_bool mz_zip_is_zip64(mz_zip_archive *pZip);
MINIZ_EXPORT size_t mz_zip_get_central_dir_size(mz_zip_archive *pZip);
MINIZ_EXPORT void *mz_zip_reader_extract_to_heap(mz_zip_archive *pZip, mz_uint file_index, size_t *pSize, mz_uint flags);
MINIZ_EXPORT void *mz_zip_reader_extract_file_to_heap(mz_zip_archive *pZip, const char *pFilename, size_t *pSize, mz_uint flags);
MINIZ_EXPORT mz_bool mz_zip_end(mz_zip_archive *pZip);
