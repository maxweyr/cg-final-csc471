#pragma once

#ifndef LAB471_SHAPE_H_INCLUDED
#define LAB471_SHAPE_H_INCLUDED

#include <string>
#include <vector>
#include <memory>
#include <map>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include "GLSL.h"

class Program;

class Shape
{
public:
    Shape();
    ~Shape();

    bool loadMesh(const std::string& path);
    void init();
    void draw(const std::shared_ptr<Program> prog) const;
    void update(float deltaTime);

    // Animation methods
    void setAnimation(int animIndex);
    int getAnimationCount() const;
    std::string getAnimationName(int index) const;
    int getCurrentAnimation() const { return currentAnimation; }

    // Transform methods
    void setPosition(const glm::vec3& pos) { position = pos; }
    void setRotation(float angle, const glm::vec3& axis) { rotationAngle = angle; rotationAxis = axis; }
    void setScale(const glm::vec3& sc) { scale = sc; }

    // Bounding box info
    glm::vec3 min = glm::vec3(0);
    glm::vec3 max = glm::vec3(0);

private:
    // Internal structures for storing mesh data
    struct Vertex {
        glm::vec3 position;
        glm::vec3 normal;
        glm::vec2 texCoords;
        glm::vec3 tangent;
        glm::vec3 bitangent;
        glm::ivec4 boneIDs;
        glm::vec4 weights;
    };

    struct Texture {
        unsigned int id;
        std::string type;
        std::string path;
    };

    struct Mesh {
        std::vector<Vertex> vertices;
        std::vector<unsigned int> indices;
        std::vector<Texture> textures;
        unsigned int VAO, VBO, EBO;
    };

    struct Bone {
        std::string name;
        int id;
        glm::mat4 offset;
    };

    struct Animation {
        std::string name;
        float duration;
        float ticksPerSecond;
        bool loop;
    };

    // Processing methods
    void processNode(aiNode* node, const aiScene* scene);
    Mesh processMesh(aiMesh* mesh, const aiScene* scene);
    std::vector<Texture> loadMaterialTextures(aiMaterial* mat, aiTextureType type, const std::string& typeName, const aiScene* scene);
    void loadAnimations(const aiScene* scene);
    void updateBoneTransforms(float time, const aiScene* scene);
    unsigned int loadTexture(const std::string& path);

    // Mesh data
    std::vector<Mesh> meshes;
    std::vector<Bone> bones;
    std::map<std::string, int> boneMap;
    std::vector<Animation> animations;
    std::vector<glm::mat4> boneTransforms;

    // Animation state
    int currentAnimation;
    float animationTime;

    // Transform state
    glm::vec3 position;
    float rotationAngle;
    glm::vec3 rotationAxis;
    glm::vec3 scale;

    // Model path info
    std::string directory;
    std::map<std::string, unsigned int> loadedTextures;

    // Assimp objects
    Assimp::Importer importer;
    const aiScene* scene;
};

#endif // LAB471_SHAPE_H_INCLUDED