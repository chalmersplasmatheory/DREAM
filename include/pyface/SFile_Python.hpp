#ifndef _DREAM_SFILE_PYTHON_H
#define _DREAM_SFILE_PYTHON_H

#ifndef PY_SSIZE_T_CLEAN
#   define PY_SSIZE_T_CLEAN
#endif
#include <Python.h>
#include <string>
#include <softlib/SFile.h>

class SFile_Python : public SFile {
private:
    PyObject *dict=nullptr;

    // Disable
    virtual void Open(const std::string&, enum sfile_mode) override;
    virtual void WriteAttribute_scalar(const std::string&, const std::string&, const double) override;
    virtual void WriteAttribute_string(const std::string&, const std::string&, const std::string&) override;

    template<typename T>
    T *getMultiArrayLinear(const std::string&, const sfilesize_t, sfilesize_t&, sfilesize_t *dims);
    template<typename T>
    T *get1DArray(const std::string&, sfilesize_t*);
    template<typename T>
    T **get2DArray(const std::string&, sfilesize_t*);

    template<typename T>
    void writeArray(const std::string&, const T*, const sfilesize_t, const sfilesize_t*, const int);

protected:
    PyObject *getDict(PyObject*, const std::string&);
    PyObject *getParentStruct(const std::string&, std::string&);
    PyObject *getObjectByName(const std::string&);

public:
    SFile_Python();
    virtual ~SFile_Python();
    
    PyObject *GetPythonDict() { return this->dict; }

    virtual bool HasVariable(const std::string&) override;

    virtual void Close() override;
    virtual void CreateStruct(const std::string&) override;
    virtual double GetAttributeScalar(const std::string&, const std::string&) override;
    virtual std::string GetAttributeString(const std::string&, const std::string&) override;
    virtual enum sfile_data_type GetDataType(const std::string&, enum sfile_data_type hint=SFILE_DATA_UNDEFINED) override;
    virtual double *GetMultiArray_linear(const std::string&, const sfilesize_t, sfilesize_t&, sfilesize_t*) override;
    virtual std::string GetString(const std::string&) override;
    virtual void WriteImage(const std::string&, const double *const*, sfilesize_t) override;
    virtual void WriteMultiArray(const std::string&, const double*, const sfilesize_t, const sfilesize_t*) override;
    virtual void WriteString(const std::string&, const std::string&) override;

    virtual double **GetDoubles(const std::string&, sfilesize_t*) override;
    virtual int32_t **GetInt32_2D(const std::string&, sfilesize_t*) override;
    virtual int64_t **GetInt64_2D(const std::string&, sfilesize_t*) override;
    virtual uint32_t **GetUInt32_2D(const std::string&, sfilesize_t*) override;
    virtual uint64_t **GetUInt64_2D(const std::string&, sfilesize_t*) override;

    virtual double *GetDoubles1D(const std::string&, sfilesize_t*) override;
    virtual int32_t *GetInt32_1D(const std::string&, sfilesize_t*) override;
    virtual int64_t *GetInt64_1D(const std::string&, sfilesize_t*) override;
    virtual uint32_t *GetUInt32_1D(const std::string&, sfilesize_t*) override;
    virtual uint64_t *GetUInt64_1D(const std::string&, sfilesize_t*) override;

    virtual void WriteArray(const std::string&, const double *const*, sfilesize_t, sfilesize_t) override;

    virtual void WriteInt32Array(const std::string&, const int32_t *const*, sfilesize_t, sfilesize_t) override;
    virtual void WriteInt64Array(const std::string&, const int64_t *const*, sfilesize_t, sfilesize_t) override;
    virtual void WriteUInt32Array(const std::string&, const uint32_t *const*, sfilesize_t, sfilesize_t) override;
    virtual void WriteUInt64Array(const std::string&, const uint64_t *const*, sfilesize_t, sfilesize_t) override;

    virtual void WriteMultiInt32Array(const std::string&, const int32_t *, const sfilesize_t, const sfilesize_t *) override;
    virtual void WriteMultiInt64Array(const std::string&, const int64_t *, const sfilesize_t, const sfilesize_t *) override;
    virtual void WriteMultiUInt32Array(const std::string&, const uint32_t *, const sfilesize_t, const sfilesize_t *) override;
    virtual void WriteMultiUInt64Array(const std::string&, const uint64_t *, const sfilesize_t, const sfilesize_t *) override;

    virtual void WriteList(const std::string&, const double*, sfilesize_t) override;

    virtual void WriteInt32List(const std::string&, const int32_t*, sfilesize_t) override;
    virtual void WriteInt64List(const std::string&, const int64_t*, sfilesize_t) override;
    virtual void WriteUInt32List(const std::string&, const uint32_t*, sfilesize_t) override;
    virtual void WriteUInt64List(const std::string&, const uint64_t*, sfilesize_t) override;
};

#endif/*_DREAM_SFILE_PYTHON_H*/
